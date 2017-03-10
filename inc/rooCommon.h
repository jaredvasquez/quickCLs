/*
 * =====================================================================================
 *
 *       Filename:  rooCommon.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  05/20/2012 11:25:13 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:
 *   Organization:
 *
 * =====================================================================================
 */

#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooAbsCategoryLValue.h"
#include "RooAbsPdf.h"
#include "RooAbsArg.h"
#include "RooCategory.h"
#include "RooCustomizer.h"
#include "RooCmdArg.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooPoisson.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRealSumPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooLinkedList.h"
#include "RooSimultaneous.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
// #include "FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooUniform.h"
#include "RooPlot.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "RooBinning.h"
#include "RooDataHist.h"

/* Put the custom class headers here */
#include "RooFormulaVarExt.h"
#include "FlexibleInterpVarExt.h"
#include "FlexibleInterpVarMkII.h"
#include "HggTwoSidedCBPdf.h"
#include "HggMG5aMCNLOLineShapePdf.h"
#include "Background.h"
