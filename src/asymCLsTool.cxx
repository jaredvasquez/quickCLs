/*
 Built into a CLs tool by Jared Vasquez (jared.vasquez@yale.edu) - March 2017
 Based on asymptotic CLs script written by Aaron Armbruster (armbrusa@umich.edu) - May 2012.
*/


#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "asymCLsTool.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

asymCLsTool::asymCLsTool() {
// ----------------------------------------------------------------------------------------------------- 

  //band configuration
  betterBands           = 1; // (recommendation = 1) improve bands by using a more appropriate asimov dataset for those points
  betterNegativeBands   = 0; // (recommendation = 0) improve also the negative bands
  profileNegativeAtZero = 0; // (recommendation = 0) profile asimov for negative bands at zero

  //other configuration
  defaultMinimizer      = "Minuit2";     // or "Minuit"
  defaultPrintLevel     = 0;             // Minuit print level
  defaultStrategy       = 0;             // Minimization strategy. 0-2. 0 = fastest, least robust. 2 = slowest, most robust
  killBelowFatal        = 1;             // In case you want to suppress RooFit warnings further, set to 1
  doBlind               = 1;             // in case your analysis is blinded
  conditionalExpected   = 1 && !doBlind; // Profiling mode for Asimov data: 0 = conditional MLEs, 1 = nominal MLEs
  doTilde               = 1;             // bound mu at zero if true and do the \tilde{q}_{mu} asymptotics
  doExp                 = 1;             // compute expected limit
  doObs                 = 1 && !doBlind; // compute observed limit
  precision             = 0.005;         // % precision in mu that defines iterative cutoff
  verbose               = 0;             // 1 = very spammy
  usePredictiveFit      = 0;             // experimental, extrapolate best fit nuisance parameters based on previous fit results
  maxRetries            = 3;             // number of minimize(fcn) retries before giving up

}


void asymCLsTool::runAsymptoticsCLs(
// ----------------------------------------------------------------------------------------------------- 
           const char* infile,
		       const char* workspaceName,
		       const char* modelConfigName,
		       const char* dataName,
		       const char* asimovDataName,
		       string folder,
		       string mass,
		       double CL,
		       TString option)
{
  TStopwatch timer;
  timer.Start();

  if (killBelowFatal) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(defaultMinimizer.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(defaultStrategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(defaultPrintLevel);
  //RooNLLVar::SetIgnoreZeroEntries(1);

  //check inputs
  TFile f(infile);
  w = (RooWorkspace*)f.Get(workspaceName);
  if (!w) {
    cout << "ERROR::Workspace: " << workspaceName << " doesn't exist!" << endl;
    return;
  }

  mc = (ModelConfig*)w->obj(modelConfigName);
  if (!mc) {
    cout << "ERROR::ModelConfig: " << modelConfigName << " doesn't exist!" << endl;
    return;
  }

  firstPOI = (RooRealVar*)mc->GetParametersOfInterest()->first();
  if(option.Contains("vbfoggf")) firstPOI=w->var("vbf_o_ggf");
  if(option.Contains("vhoggf")) firstPOI=w->var("vh_o_ggf");
  if(option.Contains("tthoggf")) firstPOI=w->var("tth_o_ggf");
  
  if(option.Contains("vbfoggf")||option.Contains("vhoggf")||option.Contains("tthoggf")){
    w->var("vbf_o_ggf")->setRange(0,100);
    w->var("vh_o_ggf")->setRange(0,100);
    w->var("tth_o_ggf")->setRange(0,100);
    w->var("mu_ggf")->setRange(0,100);
  }
  
  cout<<"POI name: "<<firstPOI->GetName()<<endl;
  data = (RooDataSet*)w->data(dataName);
  if (!data) {
    cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
    return;
  }

  mc->GetParametersOfInterest()->Print("v");

  //RooAbsPdf* pdf = mc->GetPdf();
  obs_nll = createNLL(data);//(RooNLLVar*)pdf->createNLL(*data);
  map_snapshots[obs_nll] = "nominalGlobs";
  map_data_nll[data] = obs_nll;

  w->saveSnapshot("nominalGlobs",*mc->GetGlobalObservables());
  w->saveSnapshot("nominalNuis",*mc->GetNuisanceParameters());

  global_status=0;
  RooDataSet* asimovData_0 = (RooDataSet*)w->data(asimovDataName);
  if (!asimovData_0) {
    asimovData_0 = makeAsimovData(conditionalExpected, obs_nll, 0);
    //asimovData_0 = makeAsimovData2((conditionalExpected ? obs_nll : (RooNLLVar*)NULL), 0., 0.);
  }

  int asimov0_status=global_status;
  
  asimov_0_nll = createNLL(asimovData_0);//(RooNLLVar*)pdf->createNLL(*asimovData_0);
  map_snapshots[asimov_0_nll] = "conditionalGlobs_0";
  map_data_nll[asimovData_0] = asimov_0_nll;
  setMu(0);
  map_muhat[asimov_0_nll] = 0;
  saveSnapshot(asimov_0_nll, 0);
  w->loadSnapshot("conditionalNuis_0");
  w->loadSnapshot("conditionalGlobs_0");
  map_nll_muhat[asimov_0_nll] = asimov_0_nll->getVal();

  target_CLs=1-CL;
  
  double med_limit = doExp ? getLimit(asimov_0_nll, 1.0) : 1.0;
  int med_status=global_status;

  double sigma = med_limit/sqrt(3.84); // pretty close
  double mu_up_p4_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 4), 1) + 4);
  double mu_up_p3_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 3), 1) + 3);
  double mu_up_p2_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 2), 1) + 2);
  double mu_up_p1_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 1), 1) + 1);
  double mu_up_n1_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-1), 1) - 1);
  double mu_up_n2_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-2), 1) - 2);
  double mu_up_n3_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-3), 1) - 3);
  double mu_up_n4_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-4), 1) - 4);

  double mu_up_p4 = mu_up_p4_approx;
  double mu_up_p3 = mu_up_p3_approx;
  double mu_up_p2 = mu_up_p2_approx;
  double mu_up_p1 = mu_up_p1_approx;
  double mu_up_n1 = mu_up_n1_approx;
  double mu_up_n2 = mu_up_n2_approx;
  double mu_up_n3 = mu_up_n3_approx;
  double mu_up_n4 = mu_up_n4_approx;

  firstPOI->setRange(-5*sigma, 5*sigma);
  map<int, int> N_status;
  if (betterBands && doExp) { // no better time than now to do this
    //find quantiles, starting with +2, since we should be at +1.96 right now

    double init_targetCLs = target_CLs;
    firstPOI->setRange(-5*sigma, 5*sigma);
    for (int N=2;N>=-2;N--) {
      if (N < 0 && !betterNegativeBands) continue;
      if (N == 0) continue;
      target_CLs=2*(1-ROOT::Math::gaussian_cdf(fabs(N))); // change this so findCrossing looks for sqrt(qmu95)=2
      if (N < 0) direction = -1;

      //get the acual value
      double NtimesSigma = getLimit(asimov_0_nll, N*med_limit/sqrt(3.84)); // use N * sigma(0) as an initial guess
      N_status[N] += global_status;
      sigma = NtimesSigma/N;
      cout << endl;
      cout << "Found N * sigma = " << N << " * " << sigma << endl;

      string muStr,muStrPr;
      w->loadSnapshot("conditionalGlobs_0");
      double pr_val = NtimesSigma;
      if (N < 0 && profileNegativeAtZero) pr_val = 0;
      RooDataSet* asimovData_N = makeAsimovData(1, asimov_0_nll, NtimesSigma, &muStr, &muStrPr, pr_val, 0);
      //RooDataSet* asimovData_N = makeAsimovData2(asimov_0_nll, NtimesSigma, pr_val, &muStr, &muStrPr);

      RooNLLVar* asimov_N_nll = createNLL(asimovData_N);//(RooNLLVar*)pdf->createNLL(*asimovData_N);
      map_data_nll[asimovData_N] = asimov_N_nll;
      map_snapshots[asimov_N_nll] = "conditionalGlobs"+muStrPr;
      w->loadSnapshot(map_snapshots[asimov_N_nll].c_str());
      w->loadSnapshot(("conditionalNuis"+muStrPr).c_str());
      setMu(NtimesSigma);

      double nll_val = asimov_N_nll->getVal();
      saveSnapshot(asimov_N_nll, NtimesSigma);
      map_muhat[asimov_N_nll] = NtimesSigma;
      if (N < 0 && doTilde) {
        setMu(0);
        firstPOI->setConstant(1);
        nll_val = getNLL(asimov_N_nll);
      }
      map_nll_muhat[asimov_N_nll] = nll_val;

      target_CLs = init_targetCLs;
      direction=1;
      double initial_guess = findCrossing(NtimesSigma/N, NtimesSigma/N, NtimesSigma);
      double limit = getLimit(asimov_N_nll, initial_guess);
      N_status[N] += global_status;

      if (N == 2) mu_up_p2 = limit;
      else if (N == 1) mu_up_p1 = limit;
      else if (N ==-1) mu_up_n1 = limit;
      else if (N ==-2) mu_up_n2 = limit;
      //return;
    }
    direction = 1;
    target_CLs = init_targetCLs;

  }

  w->loadSnapshot("conditionalNuis_0");
  double obs_limit = doObs ? getLimit(obs_nll, med_limit) : 0;
  int obs_status=global_status;

  bool hasFailures = false;
  if (obs_status != 0 || med_status != 0 || asimov0_status != 0) hasFailures = true;
  for (map<int, int>::iterator itr=N_status.begin();itr!=N_status.end();itr++) {
    if (itr->second != 0) hasFailures = true;
  }

  if (hasFailures) {
    cout << "--------------------------------" << endl;
    cout << "Unresolved fit failures detected" << endl;
    cout << "Asimov0:  " << asimov0_status << endl;
    for (map<int, int>::iterator itr=N_status.begin();itr!=N_status.end();itr++) {
      cout << "+" << itr->first << "sigma:  " << itr->first << endl;
    }
    cout << "Median:   " << med_status << endl;
    cout << "Observed: " << obs_status << endl;
    cout << "--------------------------------" << endl;
  }

  if (betterBands) cout << "Guess for bands" << endl;
  cout << "+2sigma:  " << mu_up_p2_approx << endl;
  cout << "+1sigma:  " << mu_up_p1_approx << endl;
  cout << "-1sigma:  " << mu_up_n1_approx << endl;
  cout << "-2sigma:  " << mu_up_n2_approx << endl;

  if (betterBands) {
    cout << endl;
    cout << "Correct bands" << endl;
    cout << "+2sigma:  " << mu_up_p2 << endl;
    cout << "+1sigma:  " << mu_up_p1 << endl;
    cout << "-1sigma:  " << mu_up_n1 << endl;
    cout << "-2sigma:  " << mu_up_n2 << endl;
  }

  cout << "Median:   " << med_limit << endl;
  cout << "Observed: " << obs_limit << endl;
  cout << endl;


  system(("mkdir -vp root/" + folder).c_str());
  system(("mkdir -vp out/" + folder).c_str());
  
  stringstream fileName;
  fileName << "root/" << folder << "/" << mass << ".root";
  TFile fout(fileName.str().c_str(),"recreate");

  ofstream ftxt(Form("out/%s/%s.txt",folder.c_str(),mass.c_str()));
  ftxt << mass
       << "\t"<<obs_limit
       << "\t"<<med_limit
       << "\t"<<mu_up_p2
       << "\t"<<mu_up_p1
       << "\t"<<mu_up_n1
       << "\t"<<mu_up_n2
       << endl;
  ftxt.close();

  TH1D* h_lim = new TH1D("limit","limit",7,0,7);
  h_lim->SetBinContent(1, obs_limit);
  h_lim->SetBinContent(2, med_limit);
  h_lim->SetBinContent(3, mu_up_p2);
  h_lim->SetBinContent(4, mu_up_p1);
  h_lim->SetBinContent(5, mu_up_n1);
  h_lim->SetBinContent(6, mu_up_n2);
  h_lim->SetBinContent(7, global_status);

  h_lim->GetXaxis()->SetBinLabel(1, "Observed");
  h_lim->GetXaxis()->SetBinLabel(2, "Expected");
  h_lim->GetXaxis()->SetBinLabel(3, "+2sigma");
  h_lim->GetXaxis()->SetBinLabel(4, "+1sigma");
  h_lim->GetXaxis()->SetBinLabel(5, "-1sigma");
  h_lim->GetXaxis()->SetBinLabel(6, "-2sigma");
  h_lim->GetXaxis()->SetBinLabel(7, "Global status"); // do something with this later

  TH1D* h_lim_old = new TH1D("limit_old","limit_old",7,0,7); // include also old approximation of bands
  h_lim_old->SetBinContent(1, obs_limit);
  h_lim_old->SetBinContent(2, med_limit);
  h_lim_old->SetBinContent(3, mu_up_p2_approx);
  h_lim_old->SetBinContent(4, mu_up_p1_approx);
  h_lim_old->SetBinContent(5, mu_up_n1_approx);
  h_lim_old->SetBinContent(6, mu_up_n2_approx);
  h_lim_old->SetBinContent(7, global_status);

  h_lim_old->GetXaxis()->SetBinLabel(1, "Observed");
  h_lim_old->GetXaxis()->SetBinLabel(2, "Expected");
  h_lim_old->GetXaxis()->SetBinLabel(3, "+2sigma");
  h_lim_old->GetXaxis()->SetBinLabel(4, "+1sigma");
  h_lim_old->GetXaxis()->SetBinLabel(5, "-1sigma");
  h_lim_old->GetXaxis()->SetBinLabel(6, "-2sigma");
  h_lim_old->GetXaxis()->SetBinLabel(7, "Global status"); 

  fout.Write();
  fout.Close();

  cout << "Finished with " << nrMinimize << " calls to minimize(nll)" << endl;
  timer.Print();
}


double asymCLsTool::getLimit(RooNLLVar* nll, double initial_guess) {
// ----------------------------------------------------------------------------------------------------- 
  cout << "------------------------" << endl;
  cout << "Getting limit for nll: " << nll->GetName() << endl;
  //get initial guess based on muhat and sigma(muhat)
  firstPOI->setConstant(0);
  global_status=0;

  if (nll == asimov_0_nll) {
    setMu(0);
    firstPOI->setConstant(1);
  }

  double muhat;
  if (map_nll_muhat.find(nll) == map_nll_muhat.end()) {
    double nll_val = getNLL(nll);
    muhat = firstPOI->getVal();
    saveSnapshot(nll, muhat);
    map_muhat[nll] = muhat;
    if (muhat < 0 && doTilde) {
      setMu(0);
      firstPOI->setConstant(1);
      nll_val = getNLL(nll);
    }

    map_nll_muhat[nll] = nll_val;
  } else {
    muhat = map_muhat[nll];
  }

  if (muhat < 0.1 || initial_guess != 0) setMu(initial_guess);
  double qmu,qmuA;
  double sigma_guess = getSigma(asimov_0_nll, firstPOI->getVal(), 0, qmu);
  double sigma_b = sigma_guess;
  double mu_guess = findCrossing(sigma_guess, sigma_b, muhat);
  double pmu = calcPmu(qmu, sigma_b, mu_guess);
  double pb = calcPb(qmu, sigma_b, mu_guess);
  double CLs = calcCLs(qmu, sigma_b, mu_guess);
  double qmu95 = getQmu95(sigma_b, mu_guess);
  setMu(mu_guess);

  cout << "Initial guess:  " << mu_guess << endl;
  cout << "Sigma(obs):     " << sigma_guess << endl;
  cout << "Sigma(mu,0):    " << sigma_b << endl;
  cout << "muhat:          " << muhat << endl;
  cout << "qmu95:          " << qmu95 << endl;
  cout << "qmu:            " << qmu << endl;
  cout << "pmu:            " << pmu << endl;
  cout << "1-pb:           " << pb << endl;
  cout << "CLs:            " << CLs << endl;
  cout << endl;

  int nrDamping = 1;
  map<double, double> guess_to_corr;
  double damping_factor = 1.0;
  //double damping_factor_pre = damping_factor;
  int nrItr = 0;
  double mu_pre = muhat;//mu_guess-10*precision*mu_guess;
  double mu_pre2 = muhat;
  while (fabs(mu_pre-mu_guess) > precision*mu_guess*direction) {
    cout << "----------------------" << endl;
    cout << "Starting iteration " << nrItr << " of " << nll->GetName() << endl;
    // do this to avoid comparing multiple minima in the conditional and unconditional fits
    if (nrItr == 0) loadSnapshot(nll, muhat);
    else if (usePredictiveFit) doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
    else loadSnapshot(asimov_0_nll, mu_pre);

    sigma_guess=getSigma(nll, mu_guess, muhat, qmu);
    saveSnapshot(nll, mu_guess);


    if (nll != asimov_0_nll) {
      if (nrItr == 0) loadSnapshot(asimov_0_nll, map_nll_muhat[asimov_0_nll]);
      else if (usePredictiveFit) {
        if (nrItr == 1) doPredictiveFit(nll, map_nll_muhat[asimov_0_nll], mu_pre, mu_guess);
        else doPredictiveFit(nll, mu_pre2, mu_pre, mu_guess);
      }
      else loadSnapshot(asimov_0_nll, mu_pre);

      sigma_b=getSigma(asimov_0_nll, mu_guess, 0, qmuA);
      saveSnapshot(asimov_0_nll, mu_guess);
    } else {
      sigma_b=sigma_guess;
      qmuA=qmu;
    }

    double corr = damping_factor*(mu_guess - findCrossing(sigma_guess, sigma_b, muhat));
    for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++) {
      if (fabs(itr->first - (mu_guess-corr)) < direction*mu_guess*0.02 && fabs(corr) > direction*mu_guess*precision) {
	      damping_factor *= 0.8;
	      cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
	      if (nrDamping++ > 10)
	      {
	        nrDamping = 1;
	        damping_factor = 1.0;
	      }
	      corr *= damping_factor;
	      break;
      }
    }

    //subtract off the difference in the new and damped correction
    guess_to_corr[mu_guess] = corr;
    mu_pre2 = mu_pre;
    mu_pre = mu_guess;
    mu_guess -= corr;

    pmu = calcPmu(qmu, sigma_b, mu_pre);
    pb = calcPb(qmu, sigma_b, mu_pre);
    CLs = calcCLs(qmu, sigma_b, mu_pre);
    qmu95 = getQmu95(sigma_b, mu_pre);

    cout << "NLL:            " << nll->GetName() << endl;
    cout << "Previous guess: " << mu_pre << endl;
    cout << "Sigma(obs):     " << sigma_guess << endl;
    cout << "Sigma(mu,0):    " << sigma_b << endl;
    cout << "muhat:          " << muhat << endl;
    cout << "pmu:            " << pmu << endl;
    cout << "1-pb:           " << pb << endl;
    cout << "CLs:            " << CLs << endl;
    cout << "qmu95:          " << qmu95 << endl;
    cout << "qmu:            " << qmu << endl;
    cout << "qmuA0:          " << qmuA << endl;
    cout << "Precision:      " << direction*mu_guess*precision << endl;
    cout << "Correction:    "  << (-corr<0?" ":"") << -corr << endl;
    cout << "New guess:      " << mu_guess << endl;
    cout << endl;

    nrItr++;
    if (nrItr > 25) {
      cout << "Infinite loop detected in getLimit(). Please intervene." << endl;
      break;
    }
  }

  cout << "Found limit for nll " << nll->GetName() << ": " << mu_guess << endl;
  cout << "Finished in " << nrItr << " iterations." << endl;
  cout << endl;
  return mu_guess;
}


double asymCLsTool::getSigma(RooNLLVar* nll, double mu, double muhat, double& qmu) {
// ----------------------------------------------------------------------------------------------------- 
  qmu = getQmu(nll, mu);
  if (verbose) cout << "qmu = " << qmu << endl;
  if (mu*direction < muhat) return fabs(mu-muhat)/sqrt(qmu);
  else if (muhat < 0 && doTilde) return sqrt(mu*mu-2*mu*muhat*direction)/sqrt(qmu);
  else return (mu-muhat)*direction/sqrt(qmu);
}


double asymCLsTool::getQmu(RooNLLVar* nll, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  double nll_muhat = map_nll_muhat[nll];
  bool isConst = firstPOI->isConstant();
  firstPOI->setConstant(1);
  setMu(mu);
  double nll_val = getNLL(nll);
  firstPOI->setConstant(isConst);
  //cout << "qmu = 2 * (" << nll_val << " - " << nll_muhat << ")" << endl;
  return 2*(nll_val-nll_muhat);
}


void asymCLsTool::saveSnapshot(RooNLLVar* nll, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  stringstream snapshotName;
  snapshotName << nll->GetName() << "_" << mu;
  w->saveSnapshot(snapshotName.str().c_str(), *mc->GetNuisanceParameters());
}


void asymCLsTool::loadSnapshot(RooNLLVar* nll, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  stringstream snapshotName;
  snapshotName << nll->GetName() << "_" << mu;
  w->loadSnapshot(snapshotName.str().c_str());
}


void asymCLsTool::doPredictiveFit(RooNLLVar* nll, double mu1, double mu2, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  if (fabs(mu2-mu) < direction*mu*precision*4) {
    loadSnapshot(nll, mu2);
    return;
  }

  // extrapolate to mu using mu1 and mu2 assuming nuis scale linear in mu
  const RooArgSet* nuis = mc->GetNuisanceParameters();
  int nrNuis = nuis->getSize();
  double* theta_mu1 = new double[nrNuis];
  double* theta_mu2 = new double[nrNuis];

  TIterator* itr = nuis->createIterator();
  RooRealVar* var=NULL;
  int counter = 0;
  loadSnapshot(nll, mu1);
  while ((var == (RooRealVar*)itr->Next())) {
    theta_mu1[counter++] = var->getVal();
  }

  itr->Reset();
  counter = 0;
  loadSnapshot(nll, mu2);
  while ((var == (RooRealVar*)itr->Next())) {
    theta_mu2[counter++] = var->getVal();
  }

  itr->Reset();
  counter = 0;
  while ((var == (RooRealVar*)itr->Next())) {
    double m = (theta_mu2[counter] - theta_mu1[counter])/(mu2-mu1);
    double b = theta_mu2[counter] - m*mu2;
    double theta_extrap = m*mu+b;
    
    var->setVal(theta_extrap);
    counter++;
  }

  delete itr;
  delete theta_mu1;
  delete theta_mu2;
}


RooNLLVar* asymCLsTool::createNLL(RooDataSet* _data) {
// ----------------------------------------------------------------------------------------------------- 
  RooArgSet nuis = *mc->GetNuisanceParameters();
  RooArgSet glob = *mc->GetGlobalObservables();
  RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*_data, Constrain(nuis), GlobalObservables(glob));
  return nll;
}


double asymCLsTool::getNLL(RooNLLVar* nll) {
// ----------------------------------------------------------------------------------------------------- 
  string snapshotName = map_snapshots[nll];
  if (snapshotName != "") w->loadSnapshot(snapshotName.c_str());
  minimize(nll);
  double val = nll->getVal();
  w->loadSnapshot("nominalGlobs");
  return val;
}


double asymCLsTool::findCrossing(double sigma_obs, double sigma, double muhat) {
// ----------------------------------------------------------------------------------------------------- 
  double mu_guess = muhat + ROOT::Math::gaussian_quantile(1-target_CLs,1)*sigma_obs*direction;
  int nrItr = 0;
  int nrDamping = 1;

  map<double, double> guess_to_corr;
  double damping_factor = 1.0;
  double mu_pre = mu_guess - 10*mu_guess*precision;
  while (fabs(mu_guess-mu_pre) > direction*mu_guess*precision) {
    double sigma_obs_extrap = sigma_obs;
    double eps = 0;
    //if (extrapolateSigma) {
    //  map<double, double> map_mu_sigma = map_nll_mu_sigma[nll];
    //}

    mu_pre = mu_guess;

    double qmu95 = getQmu95(sigma, mu_guess);
    double qmu;
    qmu = 1./sigma_obs_extrap/sigma_obs_extrap*(mu_guess-muhat)*(mu_guess-muhat);
    if (muhat < 0 && doTilde) qmu = 1./sigma_obs_extrap/sigma_obs_extrap*(mu_guess*mu_guess-2*mu_guess*muhat);

    double dqmu_dmu = 2*(mu_guess-muhat)/sigma_obs_extrap/sigma_obs_extrap - 2*qmu*eps;

    double corr = damping_factor*(qmu-qmu95)/dqmu_dmu;
    for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++) {
      if (fabs(itr->first - mu_guess) < direction*mu_guess*precision) {
	      damping_factor *= 0.8;
	      if (verbose) cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
	      if (nrDamping++ > 10) {
	        nrDamping = 1;
	        damping_factor = 1.0;
	      }
	      corr *= damping_factor;
	      break;
      }
    }
    guess_to_corr[mu_guess] = corr;

    mu_guess = mu_guess - corr;
    nrItr++;
    if (nrItr > 100) {
      cout << "Infinite loop detected in findCrossing. Please intervene." << endl;
      exit(1);
    }
    if (verbose) cout << "mu_guess = " << mu_guess << ", mu_pre = " << mu_pre << ", qmu = " << qmu << ", qmu95 = " << qmu95 << ", sigma_obs_extrap = " << sigma_obs_extrap << ", sigma = " << sigma << ", direction*mu*prec = " << direction*mu_guess*precision << endl;
  }

  return mu_guess;
}


void asymCLsTool::setMu(double mu) {
// ----------------------------------------------------------------------------------------------------- 
  if (mu != mu) {
    cout << "ERROR::POI gave nan. Please intervene." << endl;
    exit(1);
  }
  if (mu > 0 && firstPOI->getMax() < mu) firstPOI->setMax(2*mu);
  if (mu < 0 && firstPOI->getMin() > mu) firstPOI->setMin(2*mu);
  firstPOI->setVal(mu);
}


double asymCLsTool::getQmu95_brute(double sigma, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  double step_size = 0.001;
  double start = step_size;
  if (mu/sigma > 0.2) start = 0;
  for (double qmu=start;qmu<20;qmu+=step_size) {
    double CLs = calcCLs(qmu, sigma, mu);
    if (CLs < target_CLs) return qmu;
  }

  return 20;
}


double asymCLsTool::getQmu95(double sigma, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  double qmu95 = 0;
  //no sane man would venture this far down into |mu/sigma|
  double target_N = ROOT::Math::gaussian_cdf(1-target_CLs,1);
  if (fabs(mu/sigma) < 0.25*target_N) {
    qmu95 = 5.83/target_N;
  } else {
    map<double, double> guess_to_corr;
    double qmu95_guess = pow(ROOT::Math::gaussian_quantile(1-target_CLs,1),2);
    int nrItr = 0;
    int nrDamping = 1;
    double damping_factor = 1.0;
    double qmu95_pre = qmu95_guess - 10*2*qmu95_guess*precision;
    while (fabs(qmu95_guess-qmu95_pre) > 2*qmu95_guess*precision) {
      qmu95_pre = qmu95_guess;
      if (verbose) {
        cout << "qmu95_guess = " << qmu95_guess << endl;
        cout << "CLs = " << calcCLs(qmu95_guess, sigma, mu) << endl;
        cout << "Derivative = " << calcDerCLs(qmu95_guess, sigma, mu) << endl;
      }

      double corr = damping_factor*(calcCLs(qmu95_guess, sigma, mu)-target_CLs)/calcDerCLs(qmu95_guess, sigma, mu);
      for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++) {
	      if (fabs(itr->first - qmu95_guess) < 2*qmu95_guess*precision) {
	        damping_factor *= 0.8;
	        if (verbose) cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
	        if (nrDamping++ > 10) {
	          nrDamping = 1;
	          damping_factor = 1.0;
	        }
	        corr *= damping_factor;
	      }
      }

      guess_to_corr[qmu95_guess] = corr;
      qmu95_guess = qmu95_guess - corr;

      if (verbose) {
        cout << "next guess = " << qmu95_guess << endl; 
        cout << "precision = " << 2*qmu95_guess*precision << endl;
        cout << endl;
      }

      nrItr++;
      if (nrItr > 200) {
        cout << "Infinite loop detected in getQmu95. Please intervene." << endl;
        exit(1);
      }
    }
    qmu95 = qmu95_guess;
  }

  if (qmu95 != qmu95) {
    qmu95 = getQmu95_brute(sigma, mu);
  }
  if (verbose) cout << "Returning qmu95 = " << qmu95 << endl;

  return qmu95;
}


double asymCLsTool::calcCLs(double qmu_tilde, double sigma, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  double pmu = calcPmu(qmu_tilde, sigma, mu);
  double pb = calcPb(qmu_tilde, sigma, mu);
  if (verbose) {
    cout << "pmu = " << pmu << endl;
    cout << "pb = " << pb << endl;
  }
  if (pb == 1) return 0.5;
  return pmu/(1-pb);
}


double asymCLsTool::calcPmu(double qmu, double sigma, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  double pmu;
  if (qmu < mu*mu/(sigma*sigma) || !doTilde) {
    pmu = 1-ROOT::Math::gaussian_cdf(sqrt(qmu));
  } else {
    pmu = 1-ROOT::Math::gaussian_cdf((qmu+mu*mu/(sigma*sigma))/(2*fabs(mu/sigma)));
  }
  if (verbose) cout << "for pmu, qmu = " << qmu << ", sigma = " << sigma<< ", mu = " << mu << ", pmu = " << pmu << endl;
  return pmu;
}


double asymCLsTool::calcPb(double qmu, double sigma, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  if (qmu < mu*mu/(sigma*sigma) || !doTilde) {
    return 1-ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qmu));
  } else {
    return 1-ROOT::Math::gaussian_cdf((mu*mu/(sigma*sigma) - qmu)/(2*fabs(mu/sigma)));
  }
}


double asymCLsTool::calcDerCLs(double qmu, double sigma, double mu) {
// ----------------------------------------------------------------------------------------------------- 
  double dpmu_dq = 0;
  double d1mpb_dq = 0;

  if (qmu < mu*mu/(sigma*sigma)) {
    double zmu = sqrt(qmu);
    dpmu_dq = -1./(2*sqrt(qmu*2*TMath::Pi()))*exp(-zmu*zmu/2);
  } else {
    double zmu = (qmu+mu*mu/(sigma*sigma))/(2*fabs(mu/sigma));
    dpmu_dq = -1./(2*fabs(mu/sigma))*1./(sqrt(2*TMath::Pi()))*exp(-zmu*zmu/2);
  }

  if (qmu < mu*mu/(sigma*sigma)) {
    double zb = fabs(mu/sigma)-sqrt(qmu);
    d1mpb_dq = -1./sqrt(qmu*2*TMath::Pi())*exp(-zb*zb/2);
  } else {
    double zb = (mu*mu/(sigma*sigma) - qmu)/(2*fabs(mu/sigma));
    d1mpb_dq = -1./(2*fabs(mu/sigma))*1./(sqrt(2*TMath::Pi()))*exp(-zb*zb/2);
  }

  double pb = calcPb(qmu, sigma, mu);
  return dpmu_dq/(1-pb)-calcCLs(qmu, sigma, mu)/(1-pb)*d1mpb_dq;
}


int asymCLsTool::minimize(RooNLLVar* nll) {
// ----------------------------------------------------------------------------------------------------- 
  nll->enableOffsetting(true);
  nrMinimize++;
  RooAbsReal* fcn = (RooAbsReal*)nll;
  return minimize(fcn);
}


int asymCLsTool::minimize(RooAbsReal* fcn) {
// ----------------------------------------------------------------------------------------------------- 
  static int nrItr = 0;
   // cout << "Starting minimization. Using these global observables" << endl;
   // mc->GetGlobalObservables()->Print("v");

  int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
  int save_strat = strat;
  RooMinimizer minim(*fcn);
  minim.setStrategy(strat);
  minim.setPrintLevel(printLevel);
  minim.setProfile();
  minim.setEps(1);

  int status = minim.minimize( ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), 
                               ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


  // up the strategy
  if (status != 0 && status != 1 && strat < 2) {
    strat++;
    cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
    minim.setStrategy(strat);
    status = minim.minimize( ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), 
                             ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  }

  if (status != 0 && status != 1 && strat < 2) {
    strat++;
    cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
    minim.setStrategy(strat);
    status = minim.minimize( ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), 
                             ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  }

  //cout << "status is " << status << endl;

  //switch minuit version and try again
  if (status != 0 && status != 1) {
    string minType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
    string newMinType;
    if (minType == "Minuit2") newMinType = "Minuit";
    else newMinType = "Minuit2";
  
    cout << "Switching minuit type from " << minType << " to " << newMinType << endl;
  
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(newMinType.c_str());
    strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    minim.setStrategy(strat);

    status = minim.minimize( ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), 
                             ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


    if (status != 0 && status != 1 && strat < 2) {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize( ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), 
                               ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    if (status != 0 && status != 1 && strat < 2) {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize( ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), 
                               ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minType.c_str());
  }

  if (status != 0 && status != 1) {
    nrItr++;
    if (nrItr > maxRetries) {
      nrItr = 0;
      global_status++;
      cout << "WARNING::Fit failure unresolved with status " << status << endl;
      return status;
    } else {
      if (nrItr == 0) { // retry with mu=0 snapshot
        w->loadSnapshot("conditionalNuis_0");
        return minimize(fcn);
      } else if (nrItr == 1) { // retry with nominal snapshot 
        w->loadSnapshot("nominalNuis");
        return minimize(fcn);
      }
    }
  }

  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);


  if (nrItr != 0) cout << "Successful fit" << endl;
  nrItr=0;
  return status;
}


void asymCLsTool::unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter)
// ----------------------------------------------------------------------------------------------------- 
{
  if (counter > 50) {
    cout << "ERROR::Couldn't unfold constraints!" << endl;
    cout << "Initial: " << endl;
    initial.Print("v");
    cout << endl;
    cout << "Final: " << endl;
    final.Print("v");
    exit(1);
  }

  TIterator* itr = initial.createIterator();
  RooAbsPdf* pdf;
  while ((pdf = (RooAbsPdf*)itr->Next())) {
    RooArgSet nuis_tmp = nuis;
    RooArgSet constraint_set(*pdf->getAllConstraints(obs, nuis_tmp, false));
    string className(pdf->ClassName());
    if (   className != "RooGaussian" 
        && className != "RooLognormal" 
        && className != "RooGamma" 
        && className != "RooPoisson" 
        && className != "RooBifurGauss")
    {
      counter++;
      unfoldConstraints(constraint_set, final, obs, nuis, counter);
    } else {
      final.add(*pdf);
    }
  }
  delete itr;
}


RooDataSet* asymCLsTool::makeAsimovData( bool doConditional, 
// ----------------------------------------------------------------------------------------------------- 
                            RooNLLVar* conditioning_nll, 
                            double mu_val, 
                            string* mu_str, 
                            string* mu_prof_str, 
                            double mu_val_profile, 
                            bool doFit)
{
  if (mu_val_profile == -999) mu_val_profile = mu_val;


  cout << "Creating asimov data at mu = " << mu_val << ", profiling at mu = " << mu_val_profile << endl;

  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  //int strat = ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
  //int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  //ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
  //RooMinuit::SetMaxIterations(10000);
  //RooMinimizer::SetMaxFunctionCalls(10000);

  ////////////////////
  //make asimov data//
  ////////////////////
  RooAbsPdf* combPdf = mc->GetPdf();

  int _printLevel = 0;

  stringstream muStr;
  muStr << setprecision(5);
  muStr << "_" << mu_val;
  if (mu_str) *mu_str = muStr.str();

  stringstream muStrProf;
  muStrProf << setprecision(5);
  muStrProf << "_" << mu_val_profile;
  if (mu_prof_str) *mu_prof_str = muStrProf.str();

  RooRealVar* mu = (RooRealVar*)mc->GetParametersOfInterest()->first();//w->var("mu");
  cout<<Option<<endl;
  if(Option.Contains("vbfoggf")){
    mu=w->var("vbf_o_ggf");
    // doFit=false;
  }
  if(Option.Contains("vhoggf")){
    // doFit=false;
    mu=w->var("vh_o_ggf");
  }
  if(Option.Contains("tthoggf")){
    // doFit=false;
    mu=w->var("tth_o_ggf");
  }
  
  // new ratios
  // ------------------------------------------
  if (Option.Contains("muRatio_VBF")) mu=w->var("muRatio_VBF");
  if (Option.Contains("muRatio_ttH")) mu=w->var("muRatio_ttH");
  if (Option.Contains("muRatio_VH"))  mu=w->var("muRatio_VH");
  if (Option.Contains("muRatio_WH"))  mu=w->var("muRatio_WH");
  if (Option.Contains("muRatio_ZH"))  mu=w->var("muRatio_ZH");

  mu->setVal(mu_val);

  RooArgSet mc_obs = *mc->GetObservables();
  RooArgSet mc_globs = *mc->GetGlobalObservables();
  RooArgSet mc_nuis = *mc->GetNuisanceParameters();

  // pair the nuisance parameter to the global observable
  RooArgSet mc_nuis_tmp = mc_nuis;
  RooArgList nui_list("ordered_nuis");
  RooArgList glob_list("ordered_globs");
  RooArgSet constraint_set_tmp(*combPdf->getAllConstraints(mc_obs, mc_nuis_tmp, false));
  RooArgSet constraint_set;
  int counter_tmp = 0;
  unfoldConstraints(constraint_set_tmp, constraint_set, mc_obs, mc_nuis_tmp, counter_tmp);

  TIterator* cIter = constraint_set.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)cIter->Next())) {
    RooAbsPdf* pdf = (RooAbsPdf*)arg;
    if (!pdf) continue;
    //cout << "Printing pdf" << endl;
    //pdf->Print();
    //cout << "Done" << endl;
    TIterator* nIter = mc_nuis.createIterator();
    RooRealVar* thisNui = NULL;
    RooAbsArg* nui_arg;
    while ((nui_arg = (RooAbsArg*)nIter->Next())) {
      if (pdf->dependsOn(*nui_arg)) {
        thisNui = (RooRealVar*)nui_arg;
        break;
      }
    }
    delete nIter;

    //RooRealVar* thisNui = (RooRealVar*)pdf->getObservables();


    //need this incase the observable isn't fundamental. 
    //in this case, see which variable is dependent on the nuisance parameter and use that.
    RooArgSet* components = pdf->getComponents();
    //cout << "\nPrinting components" << endl;
    //components->Print();
    //cout << "Done" << endl;
    components->remove(*pdf);
    if (components->getSize()) {
      TIterator* itr1 = components->createIterator();
      RooAbsArg* arg1;
      while ((arg1 = (RooAbsArg*)itr1->Next())) {
	      TIterator* itr2 = components->createIterator();
	      RooAbsArg* arg2;
	      while ((arg2 = (RooAbsArg*)itr2->Next())) {
	        if (arg1 == arg2) continue;
	        if (arg2->dependsOn(*arg1)) {
	          components->remove(*arg1);
	        }
	      }
	      delete itr2;
      }
      delete itr1;
    }
    if (components->getSize() > 1) {
      cout << "ERROR::Couldn't isolate proper nuisance parameter" << endl;
      return NULL;
    } else if (components->getSize() == 1) {
      thisNui = (RooRealVar*)components->first();
    }



    TIterator* gIter = mc_globs.createIterator();
    RooRealVar* thisGlob = NULL;
    RooAbsArg* glob_arg;
    while ((glob_arg = (RooAbsArg*)gIter->Next())) {
      if (pdf->dependsOn(*glob_arg)) {
        thisGlob = (RooRealVar*)glob_arg;
        break;
      }
    }
    delete gIter;

    if (!thisNui || !thisGlob) {
      cout << "WARNING::Couldn't find nui or glob for constraint: " << pdf->GetName() << endl;
      //return;
      continue;
    }

    if (_printLevel >= 1) {
      cout << "Pairing nui: " << thisNui->GetName() << ", with glob: " 
          << thisGlob->GetName() << ", from constraint: " << pdf->GetName() << endl;
    }

    nui_list.add(*thisNui);
    glob_list.add(*thisGlob);

    //cout << "\nPrinting Nui/glob" << endl;
    //thisNui->Print();
    //cout << "Done nui" << endl;
    //thisGlob->Print();
    //cout << "Done glob" << endl;
  }
  delete cIter;

  //save the snapshots of nominal parameters, but only if they're not already saved
  w->saveSnapshot("tmpGlobs",*mc->GetGlobalObservables());
  w->saveSnapshot("tmpNuis",*mc->GetNuisanceParameters());
  if (!w->loadSnapshot("nominalGlobs")) {
    cout << "nominalGlobs doesn't exist. Saving snapshot." << endl;
    w->saveSnapshot("nominalGlobs",*mc->GetGlobalObservables());
  }
  else w->loadSnapshot("tmpGlobs");
  if (!w->loadSnapshot("nominalNuis")) {
    cout << "nominalNuis doesn't exist. Saving snapshot." << endl;
    w->saveSnapshot("nominalNuis",*mc->GetNuisanceParameters());
  }
  else w->loadSnapshot("tmpNuis");

  RooArgSet nuiSet_tmp(nui_list);

  mu->setVal(mu_val_profile);
  mu->setConstant(1);
  //int status = 0;
  if (doConditional && doFit) {
    mc->GetParametersOfInterest()->Print("v");
    mc->GetNuisanceParameters()->Print("v");
    mc->GetGlobalObservables()->Print("v");

    minimize(conditioning_nll);
    // cout << "Using globs for minimization" << endl;
    // mc->GetGlobalObservables()->Print("v");
    // cout << "Starting minimization.." << endl;
    // RooAbsReal* nll;
    // if (!(nll = map_data_nll[combData])) nll = combPdf->createNLL(*combData, RooFit::Constrain(nuiSet_tmp));
    // RooMinimizer minim(*nll);
    // minim.setStrategy(0);
    // minim.setPrintLevel(1);
    // status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    // if (status != 0)
    // {
    //   cout << "Fit failed for mu = " << mu->getVal() << " with status " << status << endl;
    // }
    // cout << "Done" << endl;

    //combPdf->fitTo(*combData,Hesse(false),Minos(false),PrintLevel(0),Extended(), Constrain(nuiSet_tmp));
  }
  mu->setConstant(0);
  mu->setVal(mu_val);



  // loop over the nui/glob list, grab the corresponding variable from 
  // the tmp ws, and set the glob to the value of the nui
  int nrNuis = nui_list.getSize();
  if (nrNuis != glob_list.getSize()) {
    cout << "ERROR::nui_list.getSize() != glob_list.getSize()!" << endl;
    return NULL;
  }

  for (int i=0; i<nrNuis; i++) {
    RooRealVar* nui = (RooRealVar*)nui_list.at(i);
    RooRealVar* glob = (RooRealVar*)glob_list.at(i);

    //cout << "nui: " << nui << ", glob: " << glob << endl;
    //cout << "Setting glob: " << glob->GetName() << ", which had previous val: " 
    //     << glob->getVal() << ", to conditional val: " << nui->getVal() << endl;

    glob->setVal(nui->getVal());
  }

  //save the snapshots of conditional parameters
  // cout << "Saving conditional snapshots" << endl;
  // cout << "Glob snapshot name = " << "conditionalGlobs"+muStrProf.str() << endl;
  // cout << "Nuis snapshot name = " << "conditionalNuis"+muStrProf.str() << endl;
  w->saveSnapshot(("conditionalGlobs"+muStrProf.str()).c_str(),*mc->GetGlobalObservables());
  w->saveSnapshot(("conditionalNuis" +muStrProf.str()).c_str(),*mc->GetNuisanceParameters());

  if (!doConditional)
  {
    w->loadSnapshot("nominalGlobs");
    w->loadSnapshot("nominalNuis");
  }

  if (_printLevel >= 1) cout << "Making asimov" << endl;
  //make the asimov data (snipped from Kyle)
  mu->setVal(mu_val);

  int iFrame=0;

  const char* weightName="weightVar";
  RooArgSet obsAndWeight;
  //cout << "adding obs" << endl;
  obsAndWeight.add(*mc->GetObservables());
  //cout << "adding weight" << endl;

  RooRealVar* weightVar = NULL;
  if (!(weightVar = w->var(weightName)))
  {
    w->import(*(new RooRealVar(weightName, weightName, 1,0,10000000)));
    weightVar = w->var(weightName);
  }
  //cout << "weightVar: " << weightVar << endl;
  obsAndWeight.add(*w->var(weightName));

  //cout << "defining set" << endl;
  w->defineSet("obsAndWeight",obsAndWeight);


  //////////////////////////////////////////////////////
  // MAKE ASIMOV DATA FOR OBSERVABLES
  //////////////////////////////////////////////////////

  // dummy var can just have one bin since it's a dummy
  //if(w->var("ATLAS_dummyX"))  w->var("ATLAS_dummyX")->setBins(1);

  //cout <<" check expectedData by category"<<endl;
  //RooDataSet* simData=NULL;
  RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());

  RooDataSet* asimovData;
  if (!simPdf) {
    // Get pdf associated with state from simpdf
    RooAbsPdf* pdftmp = mc->GetPdf();//simPdf->getPdf(channelCat->getLabel()) ;
	
    // Generate observables defined by the pdf associated with this state
    RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

    if (_printLevel >= 1) {
      obstmp->Print();
    }

    asimovData = new RooDataSet( ("asimovData"+muStr.str()).c_str(),
                                 ("asimovData"+muStr.str()).c_str(),
                                 RooArgSet(obsAndWeight),
                                 WeightVar(*weightVar));

    RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
    double expectedEvents = pdftmp->expectedEvents(*obstmp);
    double thisNorm = 0;
    for(int jj=0; jj<thisObs->numBins(); ++jj){
      thisObs->setBin(jj);

      thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
      if (thisNorm*expectedEvents <= 0) {
        cout << "WARNING::Detected bin with zero expected events (" << thisNorm*expectedEvents 
             << ") ! Please check your inputs. Obs = " << thisObs->GetName() << ", bin = " << jj << endl;
      }

      if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18)) {
        asimovData->add(*mc->GetObservables(), thisNorm*expectedEvents);
      }
    }
    
    if (_printLevel >= 1) {
      asimovData->Print();
      cout <<"sum entries "<<asimovData->sumEntries()<<endl;
    }

    if(asimovData->sumEntries()!=asimovData->sumEntries()){
      cout << "sum entries is nan"<<endl;
      exit(1);
    }

    //((RooRealVar*)obstmp->first())->Print();
    //cout << "expected events " << pdftmp->expectedEvents(*obstmp) << endl;
     
    w->import(*asimovData);

    if (_printLevel >= 1) {
      asimovData->Print();
      cout << endl;
    }

  } else {
    map<string, RooDataSet*> asimovDataMap;
    
    //try fix for sim pdf
    RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();//(RooCategory*)w->cat("master_channel");//(RooCategory*) (&simPdf->indexCat());
    //    TIterator* iter = simPdf->indexCat().typeIterator() ;
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;
    int nrIndices = 0;
    while((tt=(RooCatType*) iter->Next())) {
      nrIndices++;
    }

    for (int i=0;i<nrIndices;i++){
      channelCat->setIndex(i);
      iFrame++;
      // Get pdf associated with state from simpdf
      RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;
	
      // Generate observables defined by the pdf associated with this state
      RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

      if (_printLevel >= 1) {
        obstmp->Print();
        cout << "on type " << channelCat->getLabel() << " " << iFrame << endl;
      }

      RooDataSet* obsDataUnbinned = new RooDataSet( Form("combAsimovData%d",iFrame),
                                                    Form("combAsimovData%d",iFrame),
                                                    RooArgSet(obsAndWeight,*channelCat),
                                                    WeightVar(*weightVar));

      RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
      double expectedEvents = pdftmp->expectedEvents(*obstmp);
      double thisNorm = 0;
      for(int jj=0; jj<thisObs->numBins(); ++jj){
        thisObs->setBin(jj);
        thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
        if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18)) 
          obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
        //cout<<"REGTEST: thisNorm="<<thisNorm<<", expectedEvents="<<expectedEvents<<", Product="<<thisNorm*expectedEvents<<endl;
      }
    
      if (_printLevel >= 1) {
        obsDataUnbinned->Print();
        cout <<"sum entries "<<obsDataUnbinned->sumEntries()<<endl;
      }

      if(obsDataUnbinned->sumEntries()!=obsDataUnbinned->sumEntries()){
        cout << "sum entries is nan"<<endl;
        exit(1);
      }

      // ((RooRealVar*)obstmp->first())->Print();
      // cout << "pdf: " << pdftmp->GetName() << endl;
      // cout << "expected events " << pdftmp->expectedEvents(*obstmp) << endl;
      // cout << "-----" << endl;

      asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;//tempData;

      if (_printLevel >= 1) {
        cout << "channel: " << channelCat->getLabel() << ", data: ";
        obsDataUnbinned->Print();
        cout << endl;
      }
    }

    asimovData = new RooDataSet( ("asimovData"+muStr.str()).c_str(),
                                 ("asimovData"+muStr.str()).c_str(),
                                 RooArgSet(obsAndWeight,*channelCat),
                                 Index(*channelCat),
                                 Import(asimovDataMap),
                                 WeightVar(*weightVar));
    w->import(*asimovData);
  }

  //bring us back to nominal for exporting
  //w->loadSnapshot("nominalNuis");
  w->loadSnapshot("nominalGlobs");

  //ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(printLevel);

  return asimovData;
}


bool asymCLsTool::checkModel(const RooStats::ModelConfig &model, bool throwOnFail) {
// ----------------------------------------------------------------------------------------------------- 
  bool ok = true; 
  std::ostringstream errors;
  std::auto_ptr<TIterator> iter;
  RooAbsPdf *pdf = model.GetPdf(); 
  if (pdf == 0) throw std::invalid_argument("Model without Pdf");
   
  RooArgSet allowedToFloat;
  
  // Check model observables
  if (model.GetObservables() == 0) {
    ok = false; errors << "ERROR: model does not define observables.\n";
    std::cout << errors.str() << std::endl;
    if (throwOnFail) throw std::invalid_argument(errors.str()); else return false;
  } else {
    allowedToFloat.add(*model.GetObservables());
  }

  // Check model parameters of interset
  if (model.GetParametersOfInterest() == 0) {
    ok = false; errors << "ERROR: model does not define parameters of interest.\n";
  } else {
    iter.reset(model.GetParametersOfInterest()->createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
      RooRealVar *v = dynamic_cast<RooRealVar *>(a);
      if (!v) { 
        errors << "ERROR: parameter of interest " << a->GetName() << " is a " << a->ClassName() << " and not a RooRealVar\n"; 
        ok = false; continue; 
      }
      if (!pdf->dependsOn(*v)) { 
        errors << "ERROR: pdf does not depend on parameter of interest " << a->GetName() << "\n"; 
        ok = false; continue; 
      }
      allowedToFloat.add(*v);
    }
  }
  
  // Check model nuisance parameters 
  if (model.GetNuisanceParameters() != 0) {
    iter.reset(model.GetNuisanceParameters()->createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
      RooRealVar *v = dynamic_cast<RooRealVar *>(a);
      if (!v) { 
        errors << "ERROR: nuisance parameter " << a->GetName() << " is a " << a->ClassName() << " and not a RooRealVar\n"; 
        ok = false; continue; 
      }
      if (v->isConstant()) { 
        errors << "ERROR: nuisance parameter " << a->GetName() << " is constant\n"; 
        ok = false; continue; 
      }
      if (!pdf->dependsOn(*v))
      {
          errors << "WARNING: pdf does not depend on nuisance parameter, removing " << a->GetName() << "\n";
          const_cast<RooArgSet*>(model.GetNuisanceParameters())->remove(*a);
          continue;
      }
      allowedToFloat.add(*v);
    }
  }

  // check model global observables 
  if (model.GetGlobalObservables() != 0) {
    iter.reset(model.GetGlobalObservables()->createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
      RooRealVar *v = dynamic_cast<RooRealVar *>(a);
      if (!v) { ok = false; 
        errors << "ERROR: global observable " << a->GetName() << " is a " << a->ClassName() << " and not a RooRealVar\n"; continue; }
      if (!v->isConstant()) { ok = false; errors << "ERROR: global observable " << a->GetName() << " is not constant\n"; continue; }
      if (!pdf->dependsOn(*v)) { errors << "WARNING: pdf does not depend on global observable " << a->GetName() << "\n"; continue; }
    }
  }

  // check the rest of the pdf
  std::auto_ptr<RooArgSet> params(pdf->getParameters(*model.GetObservables()));
  iter.reset(params->createIterator());
  for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
    if (a->isConstant() || allowedToFloat.contains(*a)) continue;
    if (a->getAttribute("flatParam")) continue;
    errors << "WARNING: pdf parameter " << a->GetName() << " (type " << a->ClassName() << ")"
            << " is not allowed to float (it's not nuisance, poi, observable or global observable)\n";
  }
  iter.reset();
  std::cout << errors.str() << std::endl;
  if (!ok && throwOnFail) throw std::invalid_argument(errors.str());
  return ok;
}



//int main(int argc, char* argv[]){
//  if(argc<4){
//    cout<<"Usage: "<<argv[0]<<" <jobname> <infile> <mass> <option>"<<endl;
//    return 0;
//  }
//  TString wname, mname, dname, pname, aname;
//  pname=argv[1];
//  TString infile=argv[2];
//  string mass=argv[3];
//  TString option=argv[4];
//  double value=0;
//  if(argc>5) value=atof(argv[5]);
//
//  bool m_highCL=option.Contains("highCL");
//  bool m_changemass=option.Contains("changemass");
//  bool m_changect=option.Contains("changect");
//
//  Option=option;
//
//  if(m_changemass){
//    if(argc<=5) abort();
//    Hmass=value;
//  }
//
//  if(m_changect){
//    if(argc<=5) abort();
//    CT=value;
//  }
//
//  // Change these variables, if you are using different setup
//  wname="combWS";
//  mname="ModelConfig";
//  aname="asimovData_0";
//  dname="obsData";
//  dname="AsimovSB"; // Jared
//
//  system("mkdir -vp out/"+pname);
//
//  if(option.Contains("binned")) dname+="binned";
//  if(option.Contains("asimovb")) dname="AsimovB";
//  if(option.Contains("asimovsb")) dname="AsimovSB";
//  if(option.Contains("combined")){
//    wname="combined";
//    mname="ModelConfig";
//    dname="obsData";
//  }
//  double CL=0.95;
//  if(m_highCL){
//    cout<<"REGTEST: calculating 99% CL limit"<<endl;
//    CL=0.99;
//  }
//
//
//  runAsymptoticsCLs(infile.Data(),wname.Data(),mname.Data(),dname.Data(),aname.Data(),pname.Data(),mass,CL,option);
//}
