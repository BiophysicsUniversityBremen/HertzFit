#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Hertzfit procedures for GitHub
//
// Manfred Radmacher
//
// University Bremen, Germany
//
// August 29, 2023
//


strconstant FOLDER_HERTZFIT	= "root:HertzFit"
strconstant FOLDER_FV	= "root:FV"
strconstant FOLDER_FORCE	= "root:FORCE"
strconstant FOLDER_Results	= "root:Results"
 
 
	StrConstant FitRangeDataStr = "Deflection;Force;Indentation"

	constant RangeWichDataDeflection = 0
	constant RangeWichDataForce = 1
	constant RangeWichDataIndentation = 2

	constant AVERAGE	= 0
	constant APPROACH	= 1
	constant RETRACT	= 2

	StrConstant WhichCurvePopStr = "Approach;Retract;"

	constant FitModelCylinder		= 1
	constant FitModelSphere		= 2
	constant FitModelCone			= 3
	constant FitModel3Pyramid		= 4
	constant FitModel4Pyramid		= 5
	constant FitModelSphereOn3Pyramid = 6
	constant FitModelSphereOn4Pyramid = 7
	constant FitModelConeBEC		= 8
	constant FitModelPyramidBEC	= 9
	constant FitModelSphereBECGarcia	= 10
	constant FitModelSphereBECChadwick = 11

	StrConstant FitModelPopStr = "Cylinder;Sphere;Cone;3Pyramid;4Pyramid;SphereOn3Pyramid;SphereOn4Pyramid;ConeBEC;PyramidBEC;SphereBECGarcia;SphereBECChadwick;"

	constant MAXBIN_SlopeHISTO = 500


	//==========================================================================================
Function HertzFitCalc_k(FitModel, Young, Poisson, PyramidAngle, Radius, FMax) // calc. the "force constant" of the sample
	//==========================================================================================
	variable FitModel, Young, Poisson, PyramidAngle, Radius, FMax
	variable ks = nan
	switch( FitModel)
		case FitModelCone:
		case FitModelConeBEC:
			ks = sqrt(8 / Pi * tan( PyramidAngle * Pi / 180) * Young  / (1 - Poisson^2 ) * FMax )
			break
		case FitModel4Pyramid:
		case FitModelPyramidBEC:
			ks = sqrt(2 * sqrt(2) * tan( PyramidAngle * Pi / 180) * Young  / (1 - Poisson^2 ) * FMax )
			break
		case FitModel3Pyramid:
			ks = sqrt(2 * 3 * sqrt(3) / 4 * tan( PyramidAngle * Pi / 180) * Young  / (1 - Poisson^2 ) * FMax )
			break
		case FitModelSphere:
		case FitModelSphereBECGarcia:
		case FitModelSphereBECChadwick:
			ks =  (6 * Young^2  / (1 - Poisson^2 )^2 * Radius *  FMax )^(1/3)
			break
		case FitModelCylinder:
			ks = 2 * Young / ( 1 - Poisson^2 ) * Radius
			break
		case FitModelSphereOn3Pyramid:
		case FitModelSphereOn4Pyramid:
			ks = 123456789 // TODO: make the correct calculation
			break
	endswitch
	return ks
End

//====================================
Function HertzFit_InitVarsAndWaves()
	//====================================
	string SavedDataFolder = GetDataFolder(1)

	NewDataFolder/O/S	$FOLDER_HERTZFIT
	KillWaves /A/Z  
	KillVariables /A/Z
	
	variable /G ForceConstant, PoissonRatio, Delta1, Delta2, DeflectionOffset, ContactPoint, ThresholdContactPoint
	variable /G Slope
	//	variable /G checkSlopePercent
	variable/G checkTiltCorr = 0
	variable/G RemoveBadCurvesFlag = 0
	variable /G HertzPrintResults = 0
	variable /G YoungModulus, YoungApproach, YoungRetract
	variable /G PyramidAngle, TipRadius, WhichCurve, FitModel, ChiSquare, Thickness
	variable /G BrushStrength, BrushLength
	variable /G CalibrationFactor = 1
	variable/G SkipPointsAppr = 0
	variable/G SkipPointsRetr = 0
	variable/G RecalibrationVal = 1.0
	variable/G IterateFlag = 0
	variable /G k_HertzApproach, k_HertzRetract
	variable /G FMax = 0
	variable /G FitRangeWhichData = 0

	ForceConstant = 0.01
	YoungModulus = 1e4
	PoissonRatio = 0.5
	Delta1 = 5e-9
	Delta2 = 100e-9
	variable /G Delta1old = Delta1
	variable /G Delta2old = Delta2

	DeflectionOffset = 0
	ContactPoint= 0
	ThresholdContactPoint = 10e-9	// 10 nm
	FitModel = FitModel4Pyramid				// 1 = Cone, 2 = Pyramid, ...
	WhichCurve = APPROACH			// 1 = Approach, 2 = Retract, 3 = Both
	PyramidAngle = 40
	TipRadius = 40e-9
	Thickness = 1e-3

	variable nof_Points = 1008

	SetDataFolder $FOLDER_FV
	variable/G ix = 0 // FVimg: x-index
	variable/G iy = 0 // FVimg: y-index
	String /G InstrumentBrand
	String /G InstrumentName
	String /G ForceMode = "Normal"
	String /G DwellMode = "NoDwell"
	variable/G AutoToggleFlag = 1

	// Step and Dwell Params
	variable /G STEPModulationOn = 0
	variable /G nSTEPs = 0

	variable /G StepAmplitude = 0
	variable /G ModulationAmplitude = 30e-9
	variable /G ModulationFrequency = 0
	variable /G StartFreq = 1
	variable /G EndFreq = 1000
	variable /G DwellDuration = 0

	SetDataFolder SavedDataFolder
	HertzFit_InitWaves(nof_Points)
End


//###########################################################################
Function HertzFit_InitWaves(nof_Points)
	//###########################################################################
	Variable nof_Points
	PauseUpdate
	string SavedDataFolder = GetDataFolder(1)
	SetDataFolder FOLDER_HERTZFIT

	if (!exists("WaveToFit") || !exists("ZWave") || !exists("SimulatedDefl") || !exists("Force") || !exists("Indentation"))
		make/O/N=(nof_Points) WaveToFit, ZWave, SimulatedDefl, Force, Indentation
	endif
	if (!exists("DeflectionRetract") || !exists("DeflectionApproach") || !exists("LVDTretract") || !exists("LVDTapproach"))
		make/O/N=(nof_Points) DeflectionRetract, DeflectionApproach, LVDTretract, LVDTapproach
	endif
	if (!exists("SimulatedForce") || !exists("SimulatedIndentation") || !exists("SimulatedZ")  || !exists("HertzResiduals"))
		make/O/N=(nof_Points) SimulatedForce, SimulatedIndentation, SimulatedZ, HertzResiduals
	endif

	if (!exists("BluntTipContactRadius") || !exists("BluntTipIndent") )
		make/O/N=(nof_Points) BluntTipContactRadius, BluntTipIndent
	endif

	make /N=4/O w_coef

	SetScale d 0,0,"m", Indentation,SimulatedDefl,ZWave,WaveToFit
	SetScale d 0,0,"m", SimulatedIndentation,SimulatedZ
	SetScale d 0,0,"m", DeflectionRetract, DeflectionApproach, LVDTretract, LVDTapproach
	SetScale d 0,0,"N", SimulatedForce,Force, HertzResiduals

	if (!exists("W_Sigma"))
		make/O/N=6 W_Sigma
	endif

	if (!exists("DeflectionOffset"))
		variable /G DeflectionOffset
	endif

	if ( !exists("CellHeight") )
		make/O/N=(1,1) CellHeight
	endif


	SetDataFolder FOLDER_FV
	if ( !exists("FVDeflection") || !exists("FVLVDT") || !exists("FVResults") || !exists("Indexes") )
		make/O/N=(1,1,1) FVDeflection, FVLVDT, Indexes, FVResults
		
		FVResults = nan
	endif


	if ( !exists("FVAmplitude") || !exists("FVPhase") || !exists("FVResults1D") )
		make/O/N=(nof_Points) FVAmplitude, FVPhase, FVResults1D
	endif

	SetDataFolder SavedDataFolder
End



//===========================================
Function HertzFitSetActFitModelPopUpStr(FitModel)
	//===========================================
	variable FitModel
	DoWindow HertzFit
	if ( V_flag )
		DoWindow/F HertzFit
		PopupMenu ModelPopUp,mode=2,popvalue=StringFromList(FitModel-1,FitModelPopStr),value= #("\""+FitModelPopStr+"\"")
	endif
End


//===========================================
Function HertzFitSetActFitRangeDataPopUp(FitRangeWhichData)
	//===========================================
	variable FitRangeWhichData
	DoWindow HertzFit
	if ( V_flag )
		DoWindow/F HertzFit
		PopupMenu RangeDataPopUp,mode=2, popvalue=StringFromList(FitRangeWhichData,FitRangeDataStr),value= #("\""+FitRangeDataStr+"\"")
	endif
End

//======================================
Function SetVarCorrectForUnitsCalc(varNum,varStr)
	//======================================
	Variable varNum
	String varStr

	variable TheNumber
	string TheUnitStr
	sscanf varStr, "%e%s", TheNumber, TheUnitStr
	
	TheNumber = str2num(varStr)
	variable TheIndex = strsearch(varStr, " ", 0)
	
	TheUnitStr = varStr[TheIndex +1, strlen(varStr) -1 ]
	
	variable multiplier
	if ( strlen(TheUnitStr) > 1 )
		if ( strlen(TheUnitStr) > 1 )
			if ( strsearch(TheUnitStr, "f", 0 ) == 0 )
				multiplier = 1e-15
			elseif ( strsearch(TheUnitStr, "p", 0 ) == 0 )
				multiplier = 1e-12
			elseif ( strsearch(TheUnitStr, "n", 0 ) == 0 )
				multiplier = 1e-9
			elseif ( strsearch(TheUnitStr, "µ", 0 ) == 0 )
				multiplier = 1e-6
			elseif ( strsearch(TheUnitStr, "m", 0 ) == 0 )
				multiplier = 1e-3
			elseif ( strsearch(TheUnitStr, "k", 0 ) == 0 )
				multiplier = 1e3
			elseif ( strsearch(TheUnitStr, "M", 0 ) == 0 )
				multiplier = 1e6
			elseif ( strsearch(TheUnitStr, "G", 0 ) == 0 )
				multiplier = 1e9
			else
				multiplier = 1
			endif
		endif
	else
		multiplier = 1

	endif
			varNum *= multiplier

	return ( varNum )
End

//#########################################################################
Function HertzFitWindowSetVar(ctrlName,varNum,varStr,varName) : SetVariableControl
	//#########################################################################
	String ctrlName
	Variable varNum
	String varStr
	String varName

	String SavedDataFolder = GetDataFolder(1)
	SetDataFolder FOLDER_HERTZFIT

	NVar TipRadius = $(FOLDER_HERTZFIT+":TipRadius")
	varNum = SetVarCorrectForUnitsCalc(varNum,varStr)
	NVAR/Z TheVariable = $varName

//	NVAR ix = $(FOLDER_FV+":ix")
//	NVAR iy = $(FOLDER_FV+":iy")

	if (NVAR_EXISTS(TheVariable) == 1)
		TheVariable = varNum
	endif

	strswitch(varName)
		case "Delta1":
		case "Delta2":
			if ( HertzFitSanityCheckDelta1Delta2(varName) )
				HertzfitDoHertzFit(APPROACH)
//				FVanaLoop(NO_1D_LOADING, ix, iy)
				//				DoFVfit("HertzFitWindowSetVar",ctrlName)
			endif
			break
		case "ForceConstant":
		case "PoissonRatio":
		case "PyramidAngle":
		case "TipRadius":
		case "ThresholdContactPoint":
		case "SkipPointsAppr":
		case "SkipPointsRetr":
		case "RecalibrationVal":
		case "Thickness":
			//			DoFVfit("HertzFitWindowSetVar",ctrlName)
				HertzfitDoHertzFit(APPROACH)
//			FVanaLoop(NO_1D_LOADING, ix, iy)
			break
		case "YoungModulus":
		case "ContactPoint":
		case "DeflectionOffset":

			wave ZWave = $(FOLDER_HERTZFIT+":ZWave")
			wave Force = $(FOLDER_HERTZFIT+":Force")
			wave Indentation = $(FOLDER_HERTZFIT+":Indentation")
			wave WaveToFit = $(FOLDER_HERTZFIT+":WaveToFit")
			NVAR ContactPoint = $(FOLDER_HERTZFIT+":ContactPoint")

			Indentation = NAN
			Indentation = ZWave  - ContactPoint - WaveToFit
			Indentation = Indentation[p] <= 0 ? NaN : Indentation[p]

			HertzFitSimulateForceCurve()
			break
	endswitch
	SetDataFolder savedDataFolder
End


//=======================================
function HertzFitSanityCheckDelta1Delta2(varName)
	//=======================================
	String varName
	NVar Delta1 = $(FOLDER_HERTZFIT+":Delta1")
	NVar Delta2 = $(FOLDER_HERTZFIT+":Delta2")
	NVar Delta1old = $(FOLDER_HERTZFIT+":Delta1old")
	NVar Delta2old = $(FOLDER_HERTZFIT+":Delta2old")
	wave Deflection = $(FOLDER_FORCE+":Deflection")
	variable DeltaMax = WaveMax(Deflection) - WaveMin(Deflection)
	variable RetValue

	strswitch(varName)
		case "Delta1":
			if ( Delta1 >= Delta2 )
				Delta1 = Delta1old
				RetValue = 0
			else
				RetValue = 1
				Delta1old = Delta1
			endif
			break
		case "Delta2":
			if ( Delta2 <= Delta1 )
				Delta2 = Delta2old
				RetValue = 0
			else
				RetValue = 1
				Delta2old = Delta2
			endif
			break
	endswitch

	return ( RetValue )
End


//#############################################################
//		determines if force curves should be excluded from analysis
//
//		started on Nov 23 2020
//
//		by Sandra, Shruti and Manfred
//
//#############################################################

Function IsBadForceCurve(DeflectionWave, ZWave)
	wave DeflectionWave, ZWave

	String fldrSav0= GetDataFolder(1)
	SetDataFolder FOLDER_HERTZFIT

	NVAR RemoveBadCurvesFlag = $(FOLDER_HERTZFIT+":RemoveBadCurvesFlag")

	Wave W_coef

	//	print "isBadForceCurve:" , NameofWave(WaveToFit)
	SetDataFolder fldrSav0

	if ( numpnts(DeflectionWave) < 10 )
		return(1)
	else

		if ( RemoveBadCurvesFlag == 1 )
			variable p1, p2, p3, p4
			variable FinalSlope, InitialSlope
			p1 = 0
			p4 = numpnts(DeflectionWave) - 1
			p3 = 0.8 * p4
			p2 = .2 * p4

			CurveFit /Q line DeflectionWave[p1,p2] /X=ZWave
			InitialSlope = w_coef[1]
			CurveFit /Q line DeflectionWave[p3,p4]  /X=ZWave
			FinalSlope = w_coef[1]

			//			print p1, p2, p3, p4, InitialSlope, FinalSlope
			//			if ( ( InitialSlope > 0 ) && ( InitialSlope > .2 * FinalSlope ) )
			if ( ( FinalSlope < 0.01 ) || ( InitialSlope > 0.1 ) || ( InitialSlope > .1 * FinalSlope ) )
				//				print "fc removed"
				return(1)
			endif


			// add extra code here
			return(0)
		else
			return(0)
		endif
	endif

	return(0)
end

Function FuncFindPoint(TheWave, TheValue, StartFromEnd, SmallerOrLarger)
	wave TheWave
	variable TheValue
	variable StartFromEnd
	variable SmallerOrLarger

	variable nPoints, i, found, ThePoint

	found = -1
	nPoints = numpnts(TheWave)
	if ( nPoints > 0 )
	i = npoints - 1
	i = 0
	do
		if ( StartFromEnd == 1 )
			ThePoint = npoints - 1 - i
		else
			ThePoint = i
		endif
		if ( SmallerOrLarger == 1 )
			if ( Thewave[ThePoint] < TheValue )
				found = ThePoint
			endif
		else
			if ( Thewave[ThePoint] > TheValue )
				found = ThePoint
			endif	
		endif

		i += 1	
	while ( ( i < npoints ) && ( found < 0 ) )
	return (found)
	else
	return(0)
	endif
end

//=========================================
Function HowManyNumbersInWave(TheWave, p1, p2)
//=========================================
	wave TheWave
	variable p1, p2
	variable HowManyNumbers = 0

	variable i

	if ( p1 > p2 )
		i = p2
		p2 = p1
		p2 = i
	endif

	if ( p1 < 0 )
		p1 = 0
	endif

	if ( p2 >= numpnts(TheWave ) )
		p2 =  numpnts(TheWave ) - 1
	endif

	for (i = p1; i <= p2; i+= 1 )
		if ( numtype(TheWave[i]) == 0 )
			HowManyNumbers += 1
		endif
	endfor
	return(HowManyNumbers)
end

constant DisplayEachFitStep = 0	// used for debugging of HertzFit

//#############################################################
Function HertzFitDoHertzFit(WhichCurveLocal)
	//#############################################################
	Variable WhichCurveLocal // 1 = Approach , 2 = Retract

	String fldrSav0= GetDataFolder(1)
	SetDataFolder FOLDER_HERTZFIT

	NVAR WhichCurve_Panel = $(FOLDER_HERTZFIT+":WhichCurve")
	NVar ForceConstant = $(FOLDER_HERTZFIT+":ForceConstant")
	NVar PoissonRatio = $(FOLDER_HERTZFIT+":PoissonRatio")
	NVAR FitRangeWhichData = $(FOLDER_HERTZFIT+":FitRangeWhichData")

	NVar Delta1 = $(FOLDER_HERTZFIT+":Delta1")
	NVar Delta2 = $(FOLDER_HERTZFIT+":Delta2")
	NVar ContactPoint = $(FOLDER_HERTZFIT+":ContactPoint")
	NVar PyramidAngle = $(FOLDER_HERTZFIT+":PyramidAngle")
	NVar TipRadius = $(FOLDER_HERTZFIT+":TipRadius")
	NVar ThresholdContactPoint = $(FOLDER_HERTZFIT+":ThresholdContactPoint")
	NVar ChiSquare = $(FOLDER_HERTZFIT+":ChiSquare")
	NVAR checkTiltCorr = $(FOLDER_HERTZFIT+":checkTiltCorr")

	NVar Thickness = $(FOLDER_HERTZFIT+":Thickness")

	NVAR FitModel = $(FOLDER_HERTZFIT+":FitModel")
	NVAR IterateFlag = $(FOLDER_HERTZFIT+":IterateFlag")
	NVar YoungModulus = $(FOLDER_HERTZFIT+":YoungModulus")
	NVar DeflectionOffset = $(FOLDER_HERTZFIT+":DeflectionOffset")
	NVar FMax = $(FOLDER_HERTZFIT+":FMax")
	//variable YoungModulus

	//variable tn = StartRunTimeMeasurement()

	variable ContactPointOffset, ForceOffset
	variable nof_Points
	variable StartFitAtContact

	// since extract_waves changes numpoints it is essential to create wave reference afterwards
	nof_Points = HertzFitGetAnaWaves(WhichCurveLocal)

	wave ZWave = $(FOLDER_HERTZFIT+":ZWave")
	wave Force = $(FOLDER_HERTZFIT+":Force")
	wave Indentation = $(FOLDER_HERTZFIT+":Indentation")
	wave WaveToFit = $(FOLDER_HERTZFIT+":WaveToFit")

	YoungModulus = 2e3
	wave/Z W_sigma
	if ( WaveExists(W_sigma) )
		wave W_sigma
		W_sigma  = nan
	endif
	ChiSquare = nan

	PauseUpdate

	if ( IsBadForceCurve(WaveToFit, ZWave) )
		YoungModulus = nan
		ContactPoint = NAN
		return ( YoungModulus )
	endif

	// guess contact point
	// the threshold should not be larger than 5 nm, even if the force curve covers a lot of deflection range
	WaveStats /Q WaveToFit
	//		variable ContactPoint_to_guess
	variable LocalThresholdPoint = ThresholdContactPoint
	if (  LocalThresholdPoint > 0.5 * ( V_Max - V_Min) )
		LocalThresholdPoint = 0.5 * ( V_Max - V_Min)
	elseif (  LocalThresholdPoint <  0.05 * ( V_Max - V_Min) )
		LocalThresholdPoint = 0.05 * ( V_Max - V_Min)
	endif

	variable IndexOfContact = FuncFindPoint(WaveToFit, LocalThresholdPoint, 1, 1)

	if (DisplayEachFitStep == 1 )
		print " "
		print "==========================================="
		print "DoHertzFit"
		print "DoHertzfit Whichcurve:", WhichCurveLocal, "Model: ", FitModel
	endif

	if ( checkTiltCorr )
		TiltCorrection(WaveToFit,IndexOfContact,WhichCurveLocal)
		IndexOfContact = FuncFindPoint(WaveToFit, LocalThresholdPoint, 1, 1)
	endif

	if ( IndexOfContact < 0 )
		YoungModulus = NaN
		ContactPoint = NAN
		HertzFitSimulateforcecurve() // => all NaN !
		return ( NaN )
	endif

	// Calc FMax now, after possible Tilt Correction
	FMax = ForceConstant * WaveMax(WaveToFit)

	// waveToFit is the deflection data extracted and an estimate for the deflection offset is already substracted
	// there is zeros at the end of force
	// I do not understand why
	ContactPoint = Zwave[IndexOfContact]
	Force = WaveToFit * ForceConstant

	HertzFitResetdraw()
	HertzFitDrawContactPoint(ContactPoint, 0, WhichCurveLocal)
	if (DisplayEachFitStep == 1 )
		print "DoHertzfit local threshold:", LocalThresholdPoint
		print "DoHertzfit threshold index and Contact Point:", IndexOfContact, ContactPoint
	endif

	// check whether piezo moves more than a micron, otherwise height of cell probably exceeds travel range
	// of piezo
	wavestats /Q Zwave
	if ( ( V_MAX - V_Min ) < 1e-09 ) //1e-6 )
		YoungModulus = Nan
		return(YoungModulus)
	endif

	HertzFitCalculateIndentation(ContactPoint)

	// if fitrange starts at contact, first a fit is done starting at 10%, then a second fit is done starting at contactpoiint
	// this is necessary to better judge the contactpoint
	variable p1, p2
	if ( delta1 == 0 )
		StartFitAtContact = 1
		delta1 = 0.1*delta2
	else
		StartFitAtContact = 0
	endif

	switch (FitRangeWhichData)
		case RangeWichDataDeflection:
			p1 = FuncFindPoint(WaveToFit, delta1, 1, 1)
			p2 = FuncFindPoint(WaveToFit, delta2, 1, 1)
			//					HertzFitDrawFitRange(p1, p2)
			break
		case RangeWichDataForce:
			p1 = FuncFindPoint(WaveToFit, delta1 / ForceConstant, 1, 1)
			p2 = FuncFindPoint(WaveToFit, delta2 / ForceConstant, 1, 1)
			//					HertzFitDrawFitRange(p1, p2)
			break
		case RangeWichDataIndentation:
			// to be on the safe side, Hertz Guess is based on a deflection range of 25% to 75%
			WaveStats /Q WaveTofit
			variable defl1 = V_min + 0.25 * ( V_max - V_min)
			variable defl2 = V_min + 0.75 * ( V_max - V_min)
			p1 = FuncFindPoint(WaveToFit, defl1, 1, 1)
			p2 = FuncFindPoint(WaveToFit, defl2, 1, 1)

			//					p1 = FuncFindPoint(Indentation, delta1, 1, 1)
			//				p2 = FuncFindPoint(Indentation, delta2, 1, 1)
			//					HertzFitDrawFitRange(p1,p2)
			break
		default:
			p1 = FuncFindPoint(WaveToFit, delta1, 1, 1)
			p2 = FuncFindPoint(WaveToFit, delta2, 1, 1)
			break
	endswitch

	// if something is weird with p1, p2 fitrange

	if ( ( p2 < 0 ) || ( ( p2 - p1 ) < 3 ) )
		YoungModulus = nan
		//				ContactPoint = nan
		//				SaveValuesForTestLoop(p1,p2,NaN,ChiSquare) // Holger 30.3.2020
		SetDataFolder fldrSav0
		return ( YoungModulus )
	endif

	// initial guess for YoungModulus and Contact
	ContactPoint += HertzFitGuessYoungAndContact(p1, p2)

	IndexOfContact =  FuncFindPoint(ZWave, ContactPoint, 0, 0)
	HertzFitDrawContactPoint(ContactPoint, 1, WhichCurveLocal)
	HertzFitCalculateIndentation(ContactPoint)

	switch ( FitModel )
		case FitModelConeBEC:		// Chadwick Model with thickness effect for a cone:
		case FitModelPyramidBEC:		// Chadwick Model with thickness effect for a pyramid:
			YoungModulus /= HertzFitBECCorrectFactor(FitModel, PyramidAngle, TipRadius, Indentation[p2], thickness)
			break
		case FitModelSphereBECChadwick:
			YoungModulus /= HertzFitBECCorrectFactor(FitModel, PyramidAngle, TipRadius, Indentation[p2], thickness)
			break
		case FitModelSphereBECGarcia:
			YoungModulus /= HertzFitBECCorrectFactor(FitModel, PyramidAngle, TipRadius, Indentation[p2], thickness)
			break
		default:
			break
	endswitch

	if (DisplayEachFitStep == 1 )
		print "DoHertzfit Guess:", IndexOfContact, ContactPoint, YoungModulus
		print "	Fitrange p1 p2, z1, z2: ", p1, p2, ZWave[p1], ZWave[p2]

		if ( ( FitModel == FitModelConeBEC ) || ( FitModel == FitModelPyramidBEC ) || ( FitModel == FitModelSphereBECChadwick ) || ( FitModel == FitModelSphereBECGarcia ) )
			print "BEC Correction Factor: ", HertzFitBECCorrectFactor(FitModel, PyramidAngle, TipRadius, Indentation[p2], thickness), YoungModulus
		endif
		HertzFitSimulateforcecurve()
		DoUpdate /W=HertzFit
	endif

	variable		V_fitOptions=4
	variable		V_FitError = 0
	CurveFit /Q/NTHR=0 line  WaveToFit[p1,p2] /X=ZWave

	Variable TheSlope = k1
	if ( TheSlope > 0.95 ) // look wether we are on substrate ...
		YoungModulus = inf
		return(YoungModulus)
	endif

	// normal fit starts here on substrate
	// first round: contact point is approximately right
	Make/D/N=7/O W_coef, Epsilon

	W_coef[0] = YoungModulus
	W_coef[1] = PoissonRatio
	W_coef[2] = PyramidAngle
	W_coef[3] = TipRadius
	W_coef[4] = 1e-12					// Offset
	W_coef[5] = 0						// ForceOffset
	W_coef[6] = 0						// thickness for BEC, first set to zero so, that CorrectionFactor is one
	// first fit does no BEC

	String FitFunc
	switch ( FitModel )
		case FitModelCone:		// Cone:
			FitFunc = "ForcevsIndHertzCone"
			break
		case FitModel4Pyramid:		// Pyramid:
			FitFunc = "ForcevsIndHertz4Pyramid"
			break
		case FitModel3Pyramid:		// Pyramid:
			FitFunc = "ForcevsIndHertz3Pyramid"
			break
		case FitModelSphere:		// Sphere:
			FitFunc = "ForcevsIndHertzSphere"
			break

		case FitModelConeBEC:		// Chadwick Model with thickness effect for a cone:
			// Chadwick Model with thickness effect for a cone:
			// first round thickness = 0 -> no correction
			FitFunc = "ForcevsIndBECCone"
			break
		case FitModelPyramidBEC:		// Chadwick Model with thickness effect for a pyramid:
			// Chadwick Model with thickness effect for a pyramid:
			// first round thickness = 0 -> no correction
			FitFunc = "ForcevsIndBECPyramid"
			break
		case FitModelCylinder:		// Cylinder Model
			FitFunc = "ForcevsIndHertzCylinder"
			break
		case FitModelSphereOn3Pyramid:
			FitFunc = "ForcevsIndHertzSphereOn3Pyramid"
			wave /Z BluntTipContactRadius = $(FOLDER_HERTZFIT+":BluntTipContactRadius")
			if ( WaveExists(BluntTipContactRadius) == 0 )
				make /N=1000 /O BluntTipContactRadius
			endif
			BluntTipContactRadius = NAN
			CalcContactForBluntTip(TipRadius, PyramidAngle, 3, Wavemax(Indentation))
			break
		case FitModelSphereOn4Pyramid:
			FitFunc = "ForcevsIndHertzSphereOn4Pyramid"
			wave /Z BluntTipContactRadius = $(FOLDER_HERTZFIT+":BluntTipContactRadius")
			if ( WaveExists(BluntTipContactRadius) == 0)
				make /N=1000 /O BluntTipContactRadius
			endif
			BluntTipContactRadius = NAN
			CalcContactForBluntTip(TipRadius, PyramidAngle, 4, Wavemax(Indentation))
			break
		case FitModelSphereBECChadwick:
			FitFunc = "FvsIndBECSphereChadwick"
			break
		case FitModelSphereBECGarcia:
			FitFunc = "FvsIndBECSphereGarcia"
			break
		default:
			DoAlert 0,"Undefined FitModel !"
			break
	endswitch

	Epsilon = 0
	Epsilon[0] = YoungModulus / 10
	Epsilon[4] = 1e-9

	// since we have introduced epsilon, we can do everything right away
	// fitting Young Modulus and Contact Point
	// in case of BEC generate Contactpoint by the appropriate Hertz Fit and do BEC fit without adjusting ContactPoint

	variable maxNaNs = p2-p1-DimSize(W_coef,0)
	WaveStats/Q/R=[p1,p2] Indentation
	V_fitError=0

	if ( ( p2-p1 ) < 4 )							// we only have two params, so at least four data points MR 8. ASug 2020
		V_fitError=1
		YoungModulus = NaN
		SetDataFolder fldrSav0
		HertzFitSimulateforcecurve() // => all NaN !

		if (DisplayEachFitStep == 1 )
			print "DoHertzfit first fit: not enough data points"
			print "	Fitrange p1 p2, z1, z2: ", p1, p2, ZWave[p1], ZWave[p2]
		endif
		return ( NaN )
	endif


	// first fit with Hertz modell
	// first fit is done with NO BEC (Thickness set to zero
	// in case of whichData is indentaion, first fit is done in deflection
	// and then a second fit is repeated
	V_fitError=0
	wavestats /Q  /R=[p1, p2] Indentation

	if ( V_npnts < 5 )
		YoungModulus = NaN
		HertzFitSimulateforcecurve() // => all NaN !
		return ( NaN )
	endif

	FuncFit/Q/H="0111011"/NTHR=0 $FitFunc W_coef Force[p1,p2] /X=Indentation /E=Epsilon

	if ( V_fitError == 0 )
		YoungModulus = W_coef[0]
		ContactPoint += W_coef[4]
		HertzFitCalculateIndentation(ContactPoint)
		IndexOfContact =  FuncFindPoint(ZWave, ContactPoint, 0, 0)
		W_coef[4] = 1e-12
	else
		YoungModulus = NaN
		//		ContactPoint = NAN				// keep contact point from guess despite error
		HertzFitSimulateforcecurve() // => all NaN !
		return ( NaN )
	endif

	if (DisplayEachFitStep == 1 )
		print "DoHertzfit first Fit:", IndexOfContact, ContactPoint, YoungModulus
		print "			Fitrange p1 p2, z1, z2: ", p1, p2, ZWave[p1], ZWave[p2]

		//						HertzFitDisplayFitStatus(12, V_fitError, IndexOfContact, ContactPoint, YoungModulus,  W_coef, V_chisq)
		HertzFitSimulateforcecurve()
		DoUpdate /W=HertzFit
	endif

	HertzFitDrawContactPoint(ContactPoint, 3, WhichCurveLocal)

	// if FitRangeWhichData is RangeWichDataIndentation, we need to adjust p1 and p2 and redo fit
	// and do a second fit in a range defined by indentation
	if ( FitRangeWhichData == RangeWichDataIndentation )
		WaveStats /Q Indentation
		if ( delta1 < V_Min )
			p1 = V_MinRowLoc
		else
			p1 = FuncFindPoint(Indentation, delta1, 1, 1)
		endif
		if ( delta2 > V_max )
			p2 = V_MaxRowLoc
		else
			p2 = FuncFindPoint(Indentation, delta2, 1, 1)
		endif
		if (DisplayEachFitStep == 1 )
			print "DoHertzfit range adjusted for indent:", IndexOfContact, ContactPoint, YoungModulus
			print "			Fitrange p1 p2, z1, z2: ", p1, p2, ZWave[p1], ZWave[p2]
		endif

		//						HertzFitDrawFitRange(p1,p2)
		if ( ( p1 > 0 ) && ( p2 > 0 ) && ( ( p2 - p1 ) > 10 ) )
			FuncFit/Q/H="0111011"/NTHR=0 $FitFunc W_coef Force[p1,p2] /X=Indentation /E=Epsilon
		else
			YoungModulus = NaN
			//							ContactPoint = NAN
			HertzFitSimulateforcecurve() // => all NaN !
			return ( NaN )
		endif

		if (DisplayEachFitStep == 1 )
			print "DoHertzfit Repeat Fit for IndentRange:", IndexOfContact, ContactPoint, YoungModulus
		endif

		if ( V_fitError == 0 )
			YoungModulus = W_coef[0]
			ContactPoint += W_coef[4]
			HertzFitCalculateIndentation(ContactPoint)
			IndexOfContact =  FuncFindPoint(ZWave, ContactPoint, 0, 0)
			W_coef[4] = 1e-12
		else
			YoungModulus = NaN
			ContactPoint = NAN
			HertzFitSimulateforcecurve() // => all NaN !
			return ( NaN )
		endif
	endif
	// endif FitRangeWhichData is RangeWichDataIndentation

	// readjust p1 and p2 due to adjustment of force offset
	if ( StartFitAtContact == 1 )
		p1 = IndexOfContact
		delta1 = 0
	else
		switch (FitRangeWhichData)
			case RangeWichDataDeflection:
				p1 = FuncFindPoint(WaveToFit, delta1, 1, 1)
				p2 = FuncFindPoint(WaveToFit, delta2, 1, 1)
				HertzFitDrawFitRange(p1,p2)
				break
			case RangeWichDataForce:
				p1 = FuncFindPoint(WaveToFit, delta1 / ForceConstant, 1, 1)
				p2 = FuncFindPoint(WaveToFit, delta2 / ForceConstant, 1, 1)
				HertzFitDrawFitRange(p1,p2)
				break
			case RangeWichDataIndentation:
				//								wavestats /Q WaveToFit
				WaveStats /Q Indentation
				if ( delta1 < V_Min )
					p1 = V_MinRowLoc
				else
					p1 = FuncFindPoint(Indentation, delta1, 1, 1)
				endif
				if ( delta2 > V_max )
					p2 = V_MaxRowLoc
				else
					p2 = FuncFindPoint(Indentation, delta2, 1, 1)
				endif
				HertzFitDrawFitRange(p1,p2)
				break
			default:
				p1 = FuncFindPoint(WaveToFit, delta1, 1, 1)
				p2 = FuncFindPoint(WaveToFit, delta2, 1, 1)
				break
		endswitch
	endif

	if (DisplayEachFitStep == 1 )
		print "DoHertzfit range adjusted after first fit completed:", IndexOfContact, ContactPoint, YoungModulus
		print "			Fitrange p1 p2, z1, z2: ", p1, p2, ZWave[p1], ZWave[p2]
	endif

	// if iterate flag
	if ( IterateFlag == 1 )
		// start iterate flag
		// readjuste the force zero and the offset of WaveToFit value using the contact point from above
		ForceOffset = mean(Force, pnt2x(Force, IndexOfContact), pnt2x(Force, IndexOfContact - 10)  )

		WaveToFit -= ForceOffset / ForceConstant
		Force -= ForceOffset
		W_coef[4] = 1e-12
		// undo contact point adjust other wise fit will not work properly
		//
		V_fitError=0
		Indentation = ZWave  - ContactPoint - WaveToFit
		Indentation = Indentation[p] < 0 ? NaN : Indentation[p]
		W_coef[4] = 0 // 25.03.2020, Index was wrong: W_coef[3] = 0

		//						maxNaNs = p2-p1-DimSize(W_coef,cROW)
		//						WaveStats/Q/R=[p1,p2] Indentation
		//					if ( V_numNaNs < maxNaNs )
		if ( HowManyNumbersInWave(Indentation, p1, p2) > 10 )

			FuncFit/Q/H="0111011"/NTHR=0 $FitFunc W_coef Force[p1,p2] /X=Indentation /E=Epsilon
			if ( V_fitError == 0 )
				YoungModulus = W_coef[0]
				ContactPoint += W_coef[4]
				HertzFitCalculateIndentation(ContactPoint)
				IndexOfContact =  FuncFindPoint(ZWave, ContactPoint, 0, 0)
				HertzFitDrawContactPoint(ContactPoint, 4, WhichCurveLocal)
				W_coef[4] = 0 // 25.03.2020, Index was wrong: W_coef[3] = 0
			else
				YoungModulus = NaN
				ContactPoint = NAN
				HertzFitSimulateforcecurve() // => all NaN !
				return ( NaN )
			endif

			if (DisplayEachFitStep == 1 )
				print "DoHertzfit iterate fit:", IndexOfContact, ContactPoint, YoungModulus
				print "			Fitrange p1 p2, z1, z2: ", p1, p2, ZWave[p1], ZWave[p2]

				//								HertzFitDisplayFitStatus(4, V_fitError, IndexOfContact, ContactPoint, YoungModulus,  W_coef, V_chisq)
				HertzFitSimulateforcecurve()
				DoUpdate /W=HertzFit
			endif
		else
			Force = NaN
			HertzFitSimulateforcecurve() // => all NaN !
		endif

		// readjuste the force zero and the offset of WaveToFit value using the contact point from above
		ForceOffset = mean(Force, pnt2x(Force, IndexOfContact), pnt2x(Force, IndexOfContact - 10)  )
		ForceOffset = 0

		WaveToFit -= ForceOffset / ForceConstant
		Force -= ForceOffset
		W_coef[4] = 0 // 25.03.2020, Index was wrong: W_coef[3] = 0
		// undo contact point adjust other wise fit will not work properly
		//
		V_fitError=0
		Indentation = ZWave  - ContactPoint - WaveToFit
		Indentation = Indentation[p] < 0 ? NaN : Indentation[p]
		W_coef[4] = 1e-12
	endif
	// end iterate flag


	if ( HowManyNumbersInWave(Indentation, p1, p2) < 10 )
		if (DisplayEachFitStep == 1 )
			print "DoHertzfit before Second Fit: not enough data points"
		endif

		YoungModulus = NaN
		HertzFitSimulateforcecurve() // => all NaN !
		return(YoungModulus)
	endif

	//						p2 = FuncFindPoint(WaveToFit, delta2, 1, 1)
	V_fitError=0
	FuncFit/Q/H="0111011"/NTHR=0 $FitFunc W_coef Force[p1,p2] /X=Indentation /E=Epsilon
	if (DisplayEachFitStep == 1 )
		print "DoHertzfit Second Fit"
		print "E, CP: ", W_coef[0], W_coef[4], ContactPoint + W_coef[4]
	endif
	// if BEC then we will make finally a fit with the correction term for thickness
	switch ( FitModel )
		case FitModelCone:		// Cone:
		case FitModel3Pyramid:	// Pyramid:
		case FitModel4Pyramid:	// Pyramid:
		case FitModelSphere:		// Sphere:
		case FitModelCylinder:	// Cylinder:
		case FitModelSphereOn3Pyramid:
		case FitModelSphereOn4Pyramid:
			break
		case FitModelConeBEC:
		case FitModelPyramidBEC:
		case FitModelSphereBECChadwick:
		case FitModelSphereBECGarcia:
			if (numtype(thickness) == 0 )
				W_coef[6] = thickness
				W_coef[0] /= HertzFitBECCorrectFactor(FitModel, PyramidAngle, TipRadius, Indentation[p2], Thickness)
				V_fitError=0
				//				FuncFit /Q/H="0111111"/NTHR=0 $FitFunc W_coef Force[p1,p2] /X=Indentation /E=Epsilon
				FuncFit /Q/H="0111011"/NTHR=0 $FitFunc W_coef Force[p1,p2] /X=Indentation /E=Epsilon
				//				FuncFit /Q/H="0111111"/NTHR=0 $FitFunc W_coef Force[p1,p2] /X=Indentation /E=Epsilon
				if (DisplayEachFitStep == 1 )
					print "DoHertzfit BEC correction"
					print "CorrectionFactor: ", HertzFitBECCorrectFactor(FitModel, PyramidAngle, TipRadius, Indentation[p2], Thickness)
					print "E, CP: ", W_coef[0], W_coef[4], ContactPoint + W_coef[4]
				endif
			else
				YoungModulus = NaN
				HertzFitSimulateforcecurve() // => all NaN !
				return ( NaN )

			endif
			break
	endswitch

	if ( V_fitError == 0 )
		YoungModulus = W_coef[0]
		ContactPoint += W_coef[4]
		HertzFitCalculateIndentation(ContactPoint)
		IndexOfContact =  FuncFindPoint(ZWave, ContactPoint, 0, 0)
		HertzFitDrawContactPoint(ContactPoint, 5, WhichCurveLocal)
		W_coef[4] = 1e-12
	else
		YoungModulus = NaN
		ContactPoint = NAN
		HertzFitSimulateforcecurve() // => all NaN !
		return ( NaN )
	endif

	if (DisplayEachFitStep == 1 )
		print "DoHertzfit second fit (pContact, zContact, E) :", IndexOfContact, ContactPoint, YoungModulus
		print "			fitrange	p1 p2, z1, z2: ", p1, p2, ZWave[p1], ZWave[p2]
		//						HertzFitDisplayFitStatus(5, V_fitError, IndexOfContact, ContactPoint, YoungModulus,  W_coef, V_chisq)
	endif

	// if we use tilt corr we have to correct DeflectionOffset
	// Normally DeflectionOffset is calculated during GuessYoungAndContact by the first 10% of data points
	// This is needed in subsequent analysis to determine the force / deflection offset
	// e.q. in Sweep Analysis, step Analysis and PowerLaw analysis
	//
	wave DeflectionApproach = $(FOLDER_HERTZFIT+":DeflectionApproach")
	wave DeflectionRetract = $(FOLDER_HERTZFIT+":DeflectionRetract")
	variable pEnd

	if ( checkTiltCorr)
		if  ( WhichCurveLocal == 1) 	// approach
			IndexOfContact =  FuncFindPoint(ZWave, ContactPoint, 0, 0)
			pEnd = IndexOfContact - 0.1 * (numpnts(DeflectionApproach) - IndexOfContact)
			DeflectionOffset =  mean(DeflectionApproach, pnt2x(DeflectionApproach, pEnd), pnt2x(DeflectionApproach, IndexOfContact) )
		else
			IndexOfContact =  FuncFindPoint(ZWave, ContactPoint, 0, 0)
			pEnd = IndexOfContact + 0.1 * (numpnts(DeflectionApproach) - IndexOfContact)
			DeflectionOffset =  mean(DeflectionRetract, pnt2x(DeflectionRetract, IndexOfContact) , pnt2x(DeflectionRetract, pEnd))

		endif
	endif

	HertzFitSimulateForceCurve()
	ChiSquare = log(V_chisq / V_npnts)
	DoUpdate /W=HertzFit

	NVAR RemoveBadCurvesFlag = $(FOLDER_HERTZFIT+":RemoveBadCurvesFlag")

	if ( RemoveBadCurvesFlag == 1 )
		variable zMin = Zwave[0]
		variable zRange = zMin - Zwave[numpnts(ZWave) - 1]
		if ( ContactPoint < ( zMin + 0.05 * ZRange ) )
			YoungModulus = NAN
			ContactPoint = NAN
		endif
	endif

	HertzFitSimulateForceCurve()

	//	SaveValuesForTestLoop(p1,p2,TheSlope,ChiSquare) // Holger 30.3.2020

	// print "Hertzfit :", WhichCurveLocal, 	YoungModulus, ContactPoint
	//StopRunTimeMeasurementAndShow(tn,"DoHertzFit()")
	SetDataFolder fldrSav0

	// a few safety measures
	// especially indentmode can get you weird results
	if ( YoungModulus < 0 )
		YoungModulus = NAN
		ContactPoint = NAN
	endif
	if (DisplayEachFitStep == 1 )
		print "DoHertzfit final result (zContact, E) :", ContactPoint, YoungModulus

	endif

	return ( YoungModulus )
end


//======================================================================
Function TiltCorrection ( WaveToFit, IndexOfContact, WhichCurve )
	//======================================================================
	Wave WaveToFit
	Variable IndexOfContact
	Variable WhichCurve
	duplicate/o WaveToFit, LineFitWave
	variable MinIndex = DimSize(WaveToFit,0)/10

	variable FitStart, FitEnd

	if ( IndexOfContact > 0 )
		if ( WhichCurve == 1 )		//Approach
			FitStart = IndexOfContact * .1
			FitEnd = IndexOfContact * .75
			if ( FitEnd < MinIndex )
				FitEnd = MinIndex
			endif

			variable		V_fitOptions=4
			variable		V_FitError = 0
			CurveFit/Q/M=2/W=0 line, WaveToFit[FitStart, FitEnd]/D=LineFitWave
		else							//Retract
			FitStart = IndexOfContact + (numpnts(WaveToFit) - IndexOfContact ) * 0.25
			FitEnd = IndexOfContact + (numpnts(WaveToFit) - IndexOfContact ) * 0.9
			V_fitOptions=4
			V_FitError = 0
			CurveFit/Q/M=2/W=0 line, WaveToFit[0, FitStart]/D=LineFitWave

		endif

		Wave LineFitWave, W_coef
		LineFitWave = W_coef[0]+W_coef[1]*x
		WaveToFit = WaveToFit - LineFitWave

		if (DisplayEachFitStep == 1 )
			print "DoHertzfit Tiltcorrection Range (p1, p2) :", 0, FitStart
		endif
	endif
	DoUpdate
End


//==============================================
Function ForcevsIndBrush(w,D) : FitFunc
	//==============================================
	Wave w
	Variable D

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(D) = A * (  (2 * L0 / D )^(9/4) - ( D / ( 2 * L0 ) )^(3/4) )
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ D
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = A
	//CurveFitDialog/ w[1] = L0

	//	return w[0] * (  (2 * w[1] / D )^(9/4) - ( D / ( 2 * w[1] ) )^(3/4) )
	return w[0] * (  ( w[1] / ( w[1] - D  ) )^(9/4) - ( ( w[1] - D ) / ( w[1] ) )^(3/4) )
End


//======================================================================
Function BrushForceFit ( Deflection, ZWave )
	//======================================================================
	Wave Deflection, ZWave
	String fldrSav0= GetDataFolder(1)
	SetDataFolder FOLDER_HERTZFIT

	NVAR ContactPoint = $(FOLDER_HERTZFIT+":ContactPoint")

	variable IndexOfContact =  FuncFindPoint(ZWave, ContactPoint, 0, 0)

	NVAR ForceConstant = $(FOLDER_HERTZFIT+":ForceConstant")
	NVAR BrushStrength	= $(FOLDER_HERTZFIT+":BrushStrength")
	NVAR BrushLength	= $(FOLDER_HERTZFIT+":BrushLength")

	duplicate/o Deflection, ForceWaveBrush
	ForceWaveBrush = ForceConstant * Deflection
	SetScale d 0,0,"N", ForceWaveBrush

	wave SimulatedDefl = $(FOLDER_HERTZFIT+":SimulatedDefl")

	variable FitRangeUpper = IndexOfContact
	do
		FitRangeUpper += 1
	while ( ( Deflection[FitRangeUpper] > SimulatedDefl[FitRangeUpper] ) && ( FitRangeUpper < numpnts(Deflection) ) )

	variable IndexOfNewContact = FuncFindPoint( ForceWaveBrush, 0, 1, 1 )

	duplicate/o $(FOLDER_HERTZFIT+":Indentation"), IndentationBrush
	variable BrushContactPoint = ZWave[IndexOfNewContact]

	IndentationBrush = ZWave - Deflection - BrushContactPoint
	duplicate /O ForceWaveBrush, ForceWaveBrushFit
	ForceWaveBrushFit = NAN
	duplicate/o ForceWaveBrushFit, DeflectionBrushFit // For display in HertzFit window
	DeflectionBrushFit = nan

	//	Cursor /W=Graph0 A, ForceWaveBrush, IndexOfNewContact
	//	Cursor /W=Graph0 B, ForceWaveBrush, FitRangeUpper
	if ( IndexOfContact > IndexOfNewContact + 10 )
		Make/o/D/N=2 W_coef
		W_coef[0] =  ForceWaveBrush[IndexOfContact] / 3
		W_coef[1] = ( IndentationBrush[IndexOfContact] -  IndentationBrush[IndexOfNewContact] ) / .3
		Make/O/T/N=4 T_Constraints
		T_Constraints[0] = {"K0 > 1e-12","K0 < 1e-6","K1 > 1e-12","K1 < 1e-4"}

		Variable V_fitOptions=4	// Suppresses Curve Fit Window
		Variable V_fitError=0		// Suppresses error messages

		FuncFit /Q/NTHR=0 ForcevsIndBrush W_coef  ForceWaveBrush[IndexOfNewContact, FitRangeUpper] /X=IndentationBrush  /C=T_Constraints

		variable UpperIndex = FuncFindPoint( IndentationBrush, 0.9 * W_coef[1], 1, 1 )
		ForceWaveBrushFit[IndexOfNewContact, UpperIndex] = ForcevsIndBrush(W_coef,  IndentationBrush[p])
		UpperIndex = FuncFindPoint( ForceWaveBrushFit, 0.5 * ForceWaveBrush[numpnts(ForceWaveBrush)-1], 1, 1 )
		DeflectionBrushFit[IndexOfNewContact, UpperIndex] = ForcevsIndBrush(W_coef,  IndentationBrush[p])
		DeflectionBrushFit /= ForceConstant
		SetScale d 0,0,"m", ForceWaveBrush
		if ( W_coef[1] < 1e-4 )
			BrushStrength = W_coef[0]
			BrushLength   = W_coef[1]
		else
			BrushStrength = NAN
			BrushLength   = NAN
		endif
	else
		BrushStrength =nan
		BrushLength   = nan

	endif
	SetDataFolder fldrSav0

End


//#############################################################
function HertzFitCalculateIndentation(ContactPoint)
	//#############################################################
	variable ContactPoint

	wave ZWave = $(FOLDER_HERTZFIT+":ZWave")
	wave Force = $(FOLDER_HERTZFIT+":Force")
	wave Indentation = $(FOLDER_HERTZFIT+":Indentation")
	wave WaveToFit = $(FOLDER_HERTZFIT+":WaveToFit")

	Indentation = NAN
	Indentation = ZWave  - ContactPoint - WaveToFit
	//	Indentation = Indentation[p] <= 0 ? NaN : Indentation[p]
	Indentation = Indentation[p] <= 0 ? 0 : Indentation[p]
end


//#############################################################
Function HertzFitSimulateForceCurve()
	//#############################################################
	String fldrSav0= GetDataFolder(1)
	SetDataFolder FOLDER_HERTZFIT
	wave ZWave = $(FOLDER_HERTZFIT+":ZWave")
	wave SimulatedDefl = $(FOLDER_HERTZFIT+":SimulatedDefl")
	wave Force = $(FOLDER_HERTZFIT+":Force")
	wave Indentation = $(FOLDER_HERTZFIT+":Indentation")
	wave SimulatedForce = $(FOLDER_HERTZFIT+":SimulatedForce")
	wave SimulatedIndentation = $(FOLDER_HERTZFIT+":SimulatedIndentation")
	wave SimulatedZ= $(FOLDER_HERTZFIT+":SimulatedZ")
	wave HertzResiduals= $(FOLDER_HERTZFIT+":HertzResiduals")

	NVar ForceConstant = $(FOLDER_HERTZFIT+":ForceConstant")
	NVar YoungModulus = $(FOLDER_HERTZFIT+":YoungModulus")
	NVar PoissonRatio = $(FOLDER_HERTZFIT+":PoissonRatio")
	NVar Delta1 = $(FOLDER_HERTZFIT+":Delta1")
	NVar Delta2 = $(FOLDER_HERTZFIT+":Delta2")
	NVar DeflectionOffset = $(FOLDER_HERTZFIT+":DeflectionOffset")
	NVar ContactPoint = $(FOLDER_HERTZFIT+":ContactPoint")
	NVar PyramidAngle = $(FOLDER_HERTZFIT+":PyramidAngle")
	NVar TipRadius = $(FOLDER_HERTZFIT+":TipRadius")
	NVar FitModel = $(FOLDER_HERTZFIT+":FitModel")
	NVar Thickness = $(FOLDER_HERTZFIT+":Thickness")

	make /o/n=7 FitParams
	if ( ( numtype(YoungModulus) == 0 ) && ( numtype(ContactPoint) == 0 ) )
		variable IndexOfContact =  FuncFindPoint(ZWave, ContactPoint, 0, 0)

		//BinarySearchValidate(IndexOfContact)
		redimension /N=(numpnts(ZWave)) SimulatedDefl, SimulatedForce, SimulatedIndentation, SimulatedZ, HertzResiduals
		//	SimulatedIndentation = nan
		// SimulatedIndentation[IndexOfContact,] = ( p - IndexOfContact)  / ( numpnts(Indentation) - IndexOfContact) * Indentation[numpnts(Indentation)-1]
		SimulatedIndentation = Indentation

		FitParams[0]= YoungModulus
		FitParams[1]= PoissonRatio
		FitParams[2]= PyramidAngle
		FitParams[3]= TipRadius
		FitParams[4]= 0    //ContactPoint
		FitParams[5]= 0    //ForceOffset
		FitParams[6]= Thickness

		if ( numtype(YoungModulus) == 0 )
			switch( FitModel)
				case FitModelCone:
					SimulatedForce = ForcevsIndHertzCone(FitParams,SimulatedIndentation)
					break
				case FitModel3Pyramid:
					SimulatedForce = ForcevsIndHertz3Pyramid(FitParams,SimulatedIndentation)
					break
				case FitModel4Pyramid:
					SimulatedForce = ForcevsIndHertz4Pyramid(FitParams,SimulatedIndentation)
					break
				case FitModelSphere:
					SimulatedForce = ForcevsIndHertzSphere(FitParams,SimulatedIndentation)
					break
				case FitModelConeBEC:
					SimulatedForce = ForcevsIndBECCone(FitParams,SimulatedIndentation)
					break
				case FitModelPyramidBEC:
					SimulatedForce = ForcevsIndBECPyramid(FitParams,SimulatedIndentation)
					break
				case FitModelCylinder:
					SimulatedForce = ForcevsIndHertzCylinder(FitParams,SimulatedIndentation)
					break
				case FitModelSphereOn4Pyramid:
					// force recalc
					//				wave BluntTipContactRadius = $(FOLDER_HERTZFIT+":BluntTipContactRadius")
					//				BluntTipContactRadius = NAN
					SimulatedForce = ForcevsIndHertzSphereOn4Pyramid(FitParams,SimulatedIndentation)
					break
				case FitModelSphereOn3Pyramid:
					// force recalc
					//				wave BluntTipContactRadius = $(FOLDER_HERTZFIT+":BluntTipContactRadius")
					//				BluntTipContactRadius = NAN
					SimulatedForce = ForcevsIndHertzSphereOn3Pyramid(FitParams,SimulatedIndentation)
					break
				case FitModelSphereBECGarcia:
					SimulatedForce = FvsIndBECSphereGarcia(FitParams,SimulatedIndentation)
					break
				case FitModelSphereBECChadwick:
					SimulatedForce = FvsIndBECSphereChadwick(FitParams,SimulatedIndentation)
					break

			endswitch
		else
			SimulatedForce = NAN
		endif

		SimulatedDefl = SimulatedForce / ForceConstant
		SimulatedZ =  ContactPoint + (SimulatedDefl + SimulatedIndentation )
		// Hertzresiduals are plotted versus SimulatedIndentation
		//    duplicate /O HertzResiduals SimulatedIndentation

		//    HertzResiduals = Force[FindPoint(Indentation, SimulatedIndentation[p],1,1)] - SimulatedForce[p]

		HertzResiduals = Force[p] - SimulatedForce[p]
	else
		SimulatedIndentation = nan
		SimulatedDefl = nan
		SimulatedZ = nan
		SimulatedForce = NAN
	endif
	SetDataFolder fldrSav0
end


//##################################################
function HertzFitGetAnaWaves(WhichCurveLocal)
	//##################################################
	Variable WhichCurveLocal

	PauseUpdate
	Wave WaveToFit = $(FOLDER_HERTZFIT+":WaveToFit")
	Wave ZWave = $(FOLDER_HERTZFIT+":ZWave")
	Wave Force = $(FOLDER_HERTZFIT+":Force")
	Wave Indentation = $(FOLDER_HERTZFIT+":Indentation")
	Wave DeflectionRetract = $(FOLDER_HERTZFIT+":DeflectionRetract")
	Wave DeflectionApproach = $(FOLDER_HERTZFIT+":DeflectionApproach")
	Wave LVDTretract = $(FOLDER_HERTZFIT+":LVDTretract")
	Wave LVDTapproach = $(FOLDER_HERTZFIT+":LVDTapproach")

	NVAR ThresholdContactPoint = $(FOLDER_HERTZFIT+":ThresholdContactPoint")

	NVAR SkipPointsAppr = $(FOLDER_HERTZFIT+":SkipPointsAppr")
	NVAR SkipPointsRetr = $(FOLDER_HERTZFIT+":SkipPointsRetr")
	variable SkipPointsLocal

	//	if ( nof_PointsDefl == 0 || nof_PointsLVDT == 0 )
	//		return -2
	//	endif
	//
	// in the new version we have indexes in the wavenote, which help splitting the force curve
	//string TheNote = note(ZSensor)
	//string IndexList = GetNoteStr("Indexes",TheNote)

	switch(WhichCurveLocal)
		case APPROACH:
			duplicate /O DeflectionApproach, WaveToFit, Force, Indentation
			duplicate /O LVDTapproach, ZWave
			SkipPointsLocal = SkipPointsAppr
			break
		case RETRACT:
			duplicate /O DeflectionRetract, WaveToFit, Force, Indentation, ZWave
			duplicate /O LVDTretract, ZWave
			WaveToFit = DeflectionRetract[numpnts(DeflectionRetract) - p - 1]
			ZWave = LVDTretract[numpnts(LVDTretract) - p - 1]
			variable xoffset = dimoffset(WaveTofit,0) + ( numpnts(WaveToFit) - 1 ) * dimdelta(WaveToFit,0)
			SetScale/P x xoffset,-dimdelta(WaveToFit,0),GetDimLabel(WaveTofit,0,-1), ZWave, WaveToFit, Force, Indentation
			SkipPointsLocal = SkipPointsRetr
			break
			//		case 3:
			//			redimension /N=(min(numpnts(DeflectionApproach), numpnts(DeflectionRetract)))  WaveToFit, ZWave,  Force, Indentation
			//			WaveToFit = (DeflectionApproach[p] + DeflectionRetract[numpnts(DeflectionApproach) - p-1]) / 2
			//			ZWave      = (LVDTapproach[p] + LVDTretract[numpnts(DeflectionApproach) - p-1]) / 2
			//			SkipPointsLocal = ???
			//			break
	endswitch
	SetScale d 0,0,"N", Force

	if ( DimSize(WaveToFit,0) < 2 || DimSize(ZWave,0) < 2 )
		return ( 0 )
	endif

	// Delete bad data at the beginning (7.2.2013, Holger):
	if ( SkipPointsLocal < numpnts(ZWave) )
		DeletePoints 0,  SkipPointsLocal, ZWave, WaveToFit, Force, Indentation
		SetScale/P x (dimoffset(WaveTofit,0) + SkipPointsLocal * dimdelta(WaveToFit,0)),dimdelta(WaveToFit,0),GetDimLabel(WaveTofit,0,-1), ZWave, WaveToFit, Force, Indentation
	endif
	HertzFitDrawSkipPointsLine(SkipPointsLocal,WhichCurveLocal)

	//// delete NANs at the beginning:
	//if ( numtype(ZWave[0] ) != 0 || numtype(WaveToFit[0] ) != 0 )
	//	do
	//		DeletePoints 0,1, ZWave, WaveToFit, Force, Indentation
	//		SetScale/P x (dimoffset(WaveTofit,0) +  dimdelta(WaveToFit,0)),dimdelta(WaveToFit,0),GetDimLabel(WaveTofit,0,-1), ZWave, WaveToFit, Force, Indentation
	//
	//	while( (numtype(ZWave[0] ) != 0  || numtype(WaveToFit[0] ) != 0 ) && numpnts(ZWave) > 1)
	//endif

	// Calculate and substract the deflection offset
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	NVAR DeflectionOffset = $(FOLDER_HERTZFIT+":DeflectionOffset")
	if (numpnts(ZWave)>1)
		DeflectionOffset =median(WaveToFit, pnt2x(WaveToFit, 0),  pnt2x(WaveToFit,0.1*Numpnts(WaveToFit)) )
		WaveToFit = WaveToFit - DeflectionOffset
	endif
	//	print WaveToFit[0],WaveToFit[1],WaveToFit[2],WaveToFit[3],WaveToFit[4]


	//	NVAR checkTiltCorr = $(FOLDER_HERTZFIT+":checkTiltCorr")
	//	if ( checkTiltCorr ) // 07.11.2017
	//		WaveStats /Q WaveToFit
	//		variable ContactPoint_to_guess
	//		if (  0.5 * ( V_Max - V_Min) < ThresholdContactPoint )
	//			 ContactPoint_to_guess = V_Min + 0.5 * ( V_Max - V_Min) //0.05*TheMax
	//		else
	//			ContactPoint_to_guess = V_Min + ThresholdContactPoint
	//		endif
	//
	//		variable IndexOfContact = FuncFindPoint(WaveToFit, ThresholdContactPoint, 1, 1)
	//		print "IndexOfContact =", IndexOfContact
	//		TiltCorrection(WaveToFit,IndexOfContact)
	//	endif


	return ( DimSize(WaveToFit,0))
end


//#############################################################
Function HertzFitDrawFitRange( PointIndex1, PointIndex2)
	//#############################################################
	variable PointIndex1, PointIndex2

	NVar ForceConstant = $(FOLDER_HERTZFIT+":ForceConstant")
	wave ZWave = $(FOLDER_HERTZFIT+":ZWave")
	wave Force = $(FOLDER_HERTZFIT+":Force")
	wave Indentation = $(FOLDER_HERTZFIT+":Indentation")
	wave WaveToFit = $(FOLDER_HERTZFIT+":WaveToFit")

	// Draw Fit range to Graph
	DoWindow HertzFit
	if ( V_flag && DimSize(WaveToFit,0) > 0 && PointIndex1 >= 0 && PointIndex1 < Numpnts(ZWave)  && PointIndex2 >= 0 && PointIndex2 < Numpnts(ZWave) )
		PauseUpdate
		//	WaveStats /Q WaveToFit
		//	SetAxis /W=HertzFit DeflFitAxis V_Min, V_Max
		//	SetAxis /W=HertzFit ForceAxis V_Min * ForceConstant , V_Max * ForceConstant

		SetDrawLayer /W=HertzFit UserFront
		variable bottom, top
		variable value1, value2
		value1 = ZWave[PointIndex1]
		value2 = ZWave[PointIndex2]

		GetAxis/Q /W=HertzFit DeflFitAxis
		bottom = V_Min
		if ( bottom < 0 )
			bottom = 0
		endif

		top = V_Max
		SetDrawEnv /W=HertzFit  xcoord= ZAxis,ycoord= DeflFitAxis, save
		SetDrawEnv /W=HertzFit linefgc= (0,40000,0),  linethick=3, save

		DrawLine /W=HertzFit value1,bottom,value1,top
		DrawLine /W=HertzFit value2,bottom,value2,top

		GetAxis/Q /W=HertzFit ForceAxis
		bottom = V_Min
		if ( bottom < 0 )
			bottom = 0
		endif
		top = V_Max
		value1 = Indentation[PointIndex1]
		value2 = Indentation[PointIndex2]
		SetDrawEnv /W=HertzFit xcoord= IndentationAxis,ycoord= ForceAxis, save
		SetDrawEnv /W=HertzFit linefgc= (0,40000,0), save
		DrawLine /W=HertzFit value1,bottom,value1,top
		DrawLine /W=HertzFit value2,bottom,value2,top
	endif
end


//#############################################################
Function HertzFitResetDraw()
	//#############################################################
	DoWindow HertzFit
	if ( V_Flag )
		SetDrawLayer /W=HertzFit /K UserFront
	endif
end

//#############################################################
Function HertzFitDrawContactPoint(ContactPoint, WhichContactPoint, WhichCurveLocal)
	//#############################################################
	variable ContactPoint, WhichContactPoint, WhichCurveLocal

	NVAR  WhichCurve = $(FOLDER_HERTZFIT+":WhichCurve")

	variable bottom, top
	variable left, right

	if ( WhichCurve == WhichCurveLocal )
		//Draw Contact Point
		DoWindow HertzFit
		if ( V_Flag )
			GetAxis/Q /W=HertzFit DeflFitAxis
			bottom = V_Min
			top = V_Max
			GetAxis/Q /W=HertzFit ZAxis
			left = V_Min
			right = V_Max
			SetDrawLayer /W=HertzFit UserFront
			SetDrawEnv /W=HertzFit save, xcoord= ZAxis,ycoord= DeflFitAxis
			SetDrawEnv /W=HertzFit save, linethick=1

			switch(WhichContactPoint)	// numeric switch
				case 0:		// execute if case matches expression
					SetDrawEnv /W=HertzFit save, linefgc= (65535,0,0)
					break						// exit from switch
				case 1:		// execute if case matches expression
					SetDrawEnv /W=HertzFit save, linefgc= (45000,0,15000)
					break						// exit from switch
				case 2:		// execute if case matches expression
					SetDrawEnv /W=HertzFit save,linefgc= (30000,0,30000)
					break
				case 3:		// execute if case matches expression
					SetDrawEnv /W=HertzFit save, linefgc= (15000,0,45000)
					break
				case 4:		// execute if case matches expression
					SetDrawEnv /W=HertzFit save,linefgc= (0,0,65535)
					break
				case 5:		// execute if case matches expression
					SetDrawEnv /W=HertzFit save,linefgc= (0,30000,30000)
					break
				default:							// optional default expression executed
					SetDrawEnv /W=HertzFit save, linefgc= (0,0,0)
			endswitch

			DrawLine /W=HertzFit ContactPoint,bottom,ContactPoint,top
			variable ypos = (top - (WhichContactPoint/20)* (top-bottom))

			SetDrawEnv /W=HertzFit save, linethick=WhichContactPoint+1
			DrawLine /W=HertzFit ContactPoint,ypos,ContactPoint+0.1*(right-left),ypos
			DoUpdate
		endif
	endif
end


//##############################################
Function HertzFitDrawSkipPointsLine(SkipPoints,WhichCurveLocal)
	//##############################################
	Variable SkipPoints, WhichCurveLocal
	string SavedDataFolder = GetDataFolder(1)
	SetDataFolder FOLDER_FORCE
	Wave Deflection, Zsensor

	variable xvalue, yvalue
	variable LeftMax, LeftMin

	if ( SkipPoints >= 0 )
		DoWindow HertzFit
		if ( V_Flag )
			GetAxis/Q /W=HertzFit Left
			LeftMin = V_Min
			LeftMax = V_Max
			SetDrawLayer /W=HertzFit UserFront
			SetDrawEnv /W=HertzFit save, xcoord= bottom,ycoord= left
			SetDrawEnv /W=HertzFit save, linethick=2

			Wave Deflection
			if ( WhichCurveLocal == APPROACH )
				xvalue = Zsensor[SkipPoints]
				yvalue  = Deflection[SkipPoints]
				if ( ( numtype(xvalue) == 0 ) && (  numtype(yvalue) == 0 ) )
					SetDrawEnv /W=HertzFit save, linefgc= (0, 0, 65535)
					DrawLine /W=HertzFit xvalue,yvalue,xvalue,LeftMax
				endif
			elseif ( WhichCurveLocal == RETRACT )
				xvalue = Zsensor[DimSize(Zsensor,0) - 1 - SkipPoints]
				yvalue  = Deflection[DimSize(Deflection,0) - 1 - SkipPoints]
				if ( ( numtype(xvalue) == 0 ) && ( numtype(yvalue) == 0 ) )
					SetDrawEnv /W=HertzFit save, linefgc= (0, 65535, 0)
					DrawLine /W=HertzFit xvalue,yvalue,xvalue,LeftMax
				endif
			endif

		endif
	endif


	//duplicate/o Zsensor, SkipPointsLineX
	//redimension/n=2 SkipPointsLineX
	//duplicate/o Deflection, SkipPointsLineY
	//redimension/n=2 SkipPointsLineY
	//SkipPointsLineX = nan
	//SkipPointsLineY = nan
	//if ( Point > 0 )
	//	DoWindow HertzFit
	//	if ( V_Flag )
	//		Wave Deflection
	//		SkipPointsLineY[0] = ( WaveMax(Deflection) + WaveMin(Deflection) ) / 2
	//		SkipPointsLineY[1] =   WaveMin(Deflection)
	//		if ( WhichCurveLocal == APPROACH )
	//			SkipPointsLineX = Zsensor[Point]
	//		elseif ( WhichCurveLocal == RETRACT )
	//			SkipPointsLineX = Zsensor[DimSize(Zsensor,cROW) - 1 - Point]
	//		endif
	//	endif
	//endif
	SetDataFolder SavedDataFolder
End


//#####################################################################################################
Function HertzFitDisplayFitStatus(StatusNumber, FitError, IndexOfContact, ContactPoint, YoungModulus,  W_coef, Chi)
	//#####################################################################################################
	variable StatusNumber, FitError, IndexOfContact, ContactPoint, YoungModulus, chi
	wave W_coef
	printf  "DoHertzfit: Status:%2d\t FitError%3d\t  Index\t%5d\t  Contact\t%10.4e\t Young\t%10.3e\t Chi\t%8.1e\r",  StatusNumber, FitError, IndexOfContact, ContactPoint, YoungModulus, Chi
	printf "\t\tW_coef:\t\t%6.1f\t%4.2f\t%4.0f\t%8.2e\t%8.2e\t%8.2e\r"  w_coef[0], w_coef[1], w_coef[2], w_coef[3], w_coef[4],  w_coef[5]
	if ( FitError != 0 )
		print "Error in Fit: ", HertzFitTranslateError(FitError)
	endif
end


//#############################################################
Function /S HertzFitTranslateError(FitError)
	//#############################################################
	variable FitError

	string ErrorStr = " "

	if ( ( FitError & 2^0) != 0)
		ErrorStr += "Any Error \r"
	endif
	if ( ( FitError & 2^1) != 0)
		ErrorStr += "Singular Matrix \r"
	endif
	if ( ( FitError & 2^2) != 0)
		ErrorStr += "Out of Memory \r"
	endif
	if ( ( FitError & 2^3) != 0)
		ErrorStr += "Function returned NAN \r"
	endif
	if ( ( FitError & 2^4) != 0)
		ErrorStr += "Fit Function requested stop \r"
	endif
	if ( ( FitError & 2^5) != 0)
		ErrorStr += "Reentrant curve fitting \r"
	endif

	return (ErrorStr)
end


//#############################################################
Function HertzFitGuessYoungAndContact(p1, p2)
	//#############################################################
	variable p1, p2				//fit range, here search range

	NVar ForceConstant = $(FOLDER_HERTZFIT+":ForceConstant")
	NVar YoungModulus = $(FOLDER_HERTZFIT+":YoungModulus")
	NVar PoissonRatio = $(FOLDER_HERTZFIT+":PoissonRatio")
	NVar Delta1 = $(FOLDER_HERTZFIT+":Delta1")
	NVar Delta2 = $(FOLDER_HERTZFIT+":Delta2")
	NVar ContactPoint = $(FOLDER_HERTZFIT+":ContactPoint")
	NVar PyramidAngle = $(FOLDER_HERTZFIT+":PyramidAngle")
	NVar TipRadius = $(FOLDER_HERTZFIT+":TipRadius")
	NVar ChiSquare = $(FOLDER_HERTZFIT+":ChiSquare")
	NVar FitModel = $(FOLDER_HERTZFIT+":FitModel")

	wave ZWave = $(FOLDER_HERTZFIT+":ZWave")
	wave Force = $(FOLDER_HERTZFIT+":Force")
	wave Indentation = $(FOLDER_HERTZFIT+":Indentation")
	wave WaveToFit = $(FOLDER_HERTZFIT+":WaveToFit")

	variable F1, F2, i1, i2,  prefactor
	variable NewContact, NewYoung

	F1 = Force[p1]
	F2 = Force[p2]
	i1 = Indentation[p1]
	i2= Indentation[p2]

	// it can happen, if ThresholdContactP is not set nicely, that we don't have indentation data all over
	// then we shall adjust p1 and p2 such that we handle the middle range of indentation data to get something
	//
	if ( ( numtype(i1) != 0 ) || ( numtype(i2) != 0 ) )
		wavestats /Q WaveToFit
		p1 = V_minRowLoc + 0.25 * ( V_maxRowLoc - V_minRowLoc )
		p2 = V_minRowLoc + 0.75 * ( V_maxRowLoc - V_minRowLoc )
		F1 = Force[p1]
		F2 = Force[p2]
		i1 = Indentation[p1]
		i2= Indentation[p2]
	endif

	switch( FitModel)
		case FitModelCone:				// cone
		case FitModelConeBEC:
			prefactor = 1 / ( 1 - PoissonRatio^2) * 2 / Pi * tan(PyramidAngle * Pi / 180)
			NewContact = ( i1 * sqrt(F2 )- i2 * sqrt(F1 ) ) / ( sqrt(F2 ) - sqrt(F1) )
			NewYoung = F1 / prefactor * 1 / (i1- NewContact )^2
			break
		case FitModel3Pyramid:			// pyramid
			prefactor = 1 / ( 1 - PoissonRatio^2) * 1 / ( 3 * sqrt(3) / 4 ) * tan(PyramidAngle * Pi / 180)
			NewContact = ( i1 * sqrt(F2 )- i2 * sqrt(F1 ) ) / ( sqrt(F2 ) - sqrt(F1) )
			NewYoung = F1 / prefactor * 1 / (i1- NewContact )^2
			break
		case FitModel4Pyramid:			// pyramid
		case FitModelPyramidBEC:
			prefactor = 1 / ( 1 - PoissonRatio^2) * 1 / sqrt(2) * tan(PyramidAngle * Pi / 180)
			NewContact = ( i1 * sqrt(F2 )- i2 * sqrt(F1 ) ) / ( sqrt(F2 ) - sqrt(F1) )
			NewYoung = F1 / prefactor * 1 / (i1- NewContact )^2
			break
		case FitModelSphere:				// sphere
		case FitModelSphereBECGarcia:
		case FitModelSphereBECChadwick:
			prefactor = 1 / ( 1 - PoissonRatio^2) * 4/3* sqrt (TipRadius)
			NewContact = ( i1 - i2 * (F1/F2)^(2/3)  ) / ( 1- (F1/F2)^(2/3)  )
			NewYoung = F1 / prefactor  / (i1- NewContact )^(3/2)
			break
		case FitModelCylinder:			// cylinder
			prefactor = 2 / ( 1 - PoissonRatio^2) * TipRadius
			NewContact = ( i1 - i2 * (F1/F2))  / ( 1- (F1/F2)  )
			NewYoung = F1 / prefactor  / (i1- NewContact )
			break
		case FitModelSphereOn3Pyramid:
		case FitModelSphereOn4Pyramid:
			// make a pyramidal guess, but correct  Young Modulus by a fudge factor of 2
			prefactor = 1 / ( 1 - PoissonRatio^2) * 1 / sqrt(2) * tan(PyramidAngle * Pi / 180)
			NewContact = ( i1 * sqrt(F2 )- i2 * sqrt(F1 ) ) / ( sqrt(F2 ) - sqrt(F1) )
			NewYoung = 0.8 * F1 / prefactor * 1 / (i1- NewContact )^2
	endswitch
	if ( numtype(NewYoung) != 0 )
		// something went wrong, return 0 since this is only an offset
		return(0)
	else
		YoungModulus = NewYoung
		return(NewContact)
	endif
end


//####################################################################
Function HertzFitBECCorrectFactor(FitModel, OpeningAngle, TipRadius, Indentation, Thickness)
	//####################################################################
	variable FitModel, OpeningAngle, TipRadius, Indentation, Thickness

	variable chi

	if ( numtype(Thickness) == 0 )
		switch (FitModel)
			case FitModelConeBEC:
			case FitModelPyramidBEC:
				if ( Thickness > 0 )
					chi =  indentation / Thickness
					return( 1 + 1.7795 * 2 *  tan(OpeningAngle * Pi / 180) / pi^2 * chi + 16 * 1.7795^2 * (tan(OpeningAngle * Pi / 180))^2 * chi^2  )
				else
					return( 1 )
				endif
			case FitModelSphereBECGarcia:
// Garcia and Garcia, 2018, “Determination of the Elastic Moduli of a Single Cell Cultured
// on a Rigid Support by Force Microscopy.” Biophysical Journal 114: 2923–2932. 
// https://doi.org/10.1016/j.bpj.2018.05.012.
				if ( Thickness > 0 )
					chi = sqrt( TipRadius * indentation ) / thickness
					if ( indentation >= thickness )
						return(NAN)
					else
						return( 1 + 1.133 * chi+ 1.497 * chi^2 + 1.469 * chi^3 + 0.755 * chi^4  )
					endif
				else
					return( 1 )
				endif
				break
			case FitModelSphereBECChadwick:
				chi = sqrt( TipRadius * indentation ) / thickness
				if ( Thickness > 0 )
					return( 1 + 1.7795 * 2 *  tan(OpeningAngle * Pi / 180) / pi^2 * chi + 16 * 1.7795^2 * (tan(OpeningAngle * Pi / 180))^2 * chi^2  )
				else
					return( 1 )
				endif
				break
			default:
				return(1)
		endswitch
	else
		return(1)
	endif

end


//#############################################################
Function ForcevsIndBECCone(w,Indentation) : FitFunc
	//#############################################################
	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable OpeningAngle = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]
	variable Thickness = w[6]

	variable CorrectionFactor

	if ( Thickness > 0 )
		//	CorrectionFactor =   1 + 1.7795 * 2 *  tan(w[2] * Pi / 180) / pi^2 * indentation / w[5]+ 16 * 1.7795^2 * (tan(w[2] * Pi / 180))^2 * ( indentation  )^2 / w[5]^2
		//    CorrectionFactor =   1 + 1.7795 * 2 *  tan(alpha * Pi / 180) / pi^2 *  ( indentation - ContactPoint ) / Thickness + 16 * 1.7795^2 * (tan(alpha * Pi / 180))^2 * (  ( indentation - ContactPoint )  )^2 / Thickness^2
		CorrectionFactor = HertzFitBECCorrectFactor(FitModelConeBEC, OpeningAngle, TipRadius, Indentation, Thickness)
	else
		CorrectionFactor = 1
	endif

	if ( ( indentation - ContactPoint ) < 0 )
		return(0)
	else
		return ( E/ ( 1 - v^2) * 2 / Pi * tan(OpeningAngle * Pi / 180) * ( indentation - ContactPoint )^2  * CorrectionFactor -  ForceOffset )
	endif
End



//#############################################################
Function ForcevsIndBECPyramid(w,Indentation) : FitFunc
	//#############################################################
	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable OpeningAngle = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]
	variable Thickness = w[6]

	variable CorrectionFactor, Result

	if ( Thickness > 0 )
		//	CorrectionFactor =   1 + 1.7795 * 2 *  tan(w[2] * Pi / 180) / pi^2 * indentation / w[5]+ 16 * 1.7795^2 * (tan(w[2] * Pi / 180))^2 * ( indentation  )^2 / w[5]^2
		// CorrectionFactor =   1 + 1.7795 * 2 *  tan(alpha * Pi / 180) / pi^2 * ( indentation - ContactPoint ) / Thickness + 16 * 1.7795^2 * (tan(alpha * Pi / 180))^2 * ( ( indentation - ContactPoint )  )^2 / Thickness^2
		CorrectionFactor = HertzFitBECCorrectFactor(FitModelPyramidBEC, OpeningAngle, TipRadius, Indentation, Thickness)
	else
		CorrectionFactor = 1
	endif

	if ( ( indentation - ContactPoint ) < 0 )
		return(0)
	else
		return( ( E / ( 1 - v^2) * 1 / sqrt(2) * tan(OpeningAngle * Pi / 180) * ( indentation - ContactPoint )^2  * CorrectionFactor -  ForceOffset ) )
	endif
End


//#############################################################
Function ForcevsIndHertzCone(w,Indentation) : FitFunc
	//#############################################################
	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable alpha = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]

	if ( ( indentation - ContactPoint ) < 0 )
		return(0)
	else
		return(  ( E / ( 1 - v^2) * 2 / Pi * tan(alpha * Pi / 180) * ( indentation - ContactPoint )^2 ) - ForceOffset)
	endif
End


//#############################################################
Function ForcevsIndHertz3Pyramid(w,Indentation) : FitFunc
	//#############################################################
	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable alpha = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]

	if ( ( indentation - ContactPoint ) < 0 )
		return(0)
	else
		return( ( E / ( 1 - v^2) * 1 / ( 3 * sqrt(3) / 4 ) * tan(alpha * Pi / 180) * ( indentation - ContactPoint )^2 ) - ForceOffset )
	endif
End



//#############################################################
Function ForcevsIndHertz4Pyramid(w,Indentation) : FitFunc
	//#############################################################
	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable alpha = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]

	if ( ( indentation - ContactPoint ) < 0 )
		return(0)
	else
		return( ( E / ( 1 - v^2) * 1 / sqrt(2) * tan(alpha * Pi / 180) * ( indentation - ContactPoint )^2 ) - ForceOffset )
	endif
End

//#############################################################
Function ForcevsIndHertzSphere(w,Indentation) : FitFunc
	//#############################################################
	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable alpha = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]

	if ( ( indentation - ContactPoint ) < 0 )
		return(0)
	else
		return E / ( 1 - v^2) * 4 / 3 * sqrt(TipRadius) * ( indentation - ContactPoint )^(3/2)  - ForceOffset
	endif

End


//#############################################################
Function FvsIndBECSphereGarcia(w,Indentation) : FitFunc
	//#############################################################
	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable OpeningAngle = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]
	variable Thickness = w[6]

	variable CorrectionFactor
	if ( Thickness > 0 )
		CorrectionFactor = HertzFitBECCorrectFactor(FitModelSphereBECGarcia, OpeningAngle, TipRadius,  indentation - ContactPoint, Thickness)
	else
		CorrectionFactor = 1
	endif

	if ( ( indentation - ContactPoint ) < 0 )
		return(0)
	else
		return E / ( 1 - v^2) * 4 / 3 * sqrt(TipRadius) * ( indentation - ContactPoint )^(3/2) * CorrectionFactor   - ForceOffset
	endif

End

//#############################################################
Function FvsIndBECSphereChadwick(w,Indentation) : FitFunc
	//#############################################################
	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable OpeningAngle = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]
	variable Thickness = w[6]

	variable CorrectionFactor
	if ( Thickness > 0 )
		CorrectionFactor = HertzFitBECCorrectFactor(FitModelSphereBECChadwick, OpeningAngle, TipRadius, Indentation, Thickness)
	else
		CorrectionFactor = 1
	endif

	if ( ( indentation - ContactPoint ) < 0 )
		return(0)
	else
		return E / ( 1 - v^2) * 4 / 3 * sqrt(TipRadius) * ( indentation - ContactPoint )^(3/2)   * CorrectionFactor   - ForceOffset
	endif

End


//#############################################################
Function ForcevsIndHertzCylinder(w,Indentation) : FitFunc
	//#############################################################
	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable alpha = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]

	if ( ( indentation - ContactPoint ) < 0 )
		return(0)
	else
		return( 2 * E / ( 1 - v^2 ) * TipRadius * ( Indentation - ContactPoint ) - ForceOffset )
	endif
End

//#################################################################
Function ForcevsIndHertzSphereOn4Pyramid(w,Indentation) : FitFunc
	//#################################################################
	//
	// this functions follows the paper by Felix Rico et al. Blunted Pyramidal Tips
	// PhysRev E 2005
	//
	// this paper does not give an analytic solution for F vs indent
	// instead the indentation is calculated from the contact radius
	// the force then needs the indnetation and the contact radius
	//
	// to speed up the process the indentation vs contact radius, which does not depend on E is
	// calculated before (in HertzFitDoHertzFit) and only used here as a LookUp Table
	//
	// The Rico model is usable for 3-sided and 4-sided pyramids -> nPyramid
	//
	// the transition point from sphere to pyramid is called b (distance from very end to plane of transition)
	// we assume here a smooth transition, i.e. the transition occurs were the angle of the sphere becomes identical
	// with the pyramid angle
	//
	// Manfred Radmacher Nov. 26, 2019
	//

	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable alpha = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]

	variable nPyramid = 4
	variable F

	wave BluntTipContactRadius = $(FOLDER_HERTZFIT+":BluntTipContactRadius")
	wave  BluntTipIndent = $(FOLDER_HERTZFIT+":BluntTipIndent")
	NVAR TransitionRadius = $(FOLDER_HERTZFIT+":TransistionRadius")

	if (numtype(Indentation) == 0 )

		//		variable TransitionRadius = CalcContactForBluntTip(TipRadius, alpha, nPyramid)

		//		variable TheContactRadius = interp(Indentation, BluntTipIndent, BluntTipContactRadius )
		variable TheContactRadius = interp(( Indentation - ContactPoint ), BluntTipIndent, BluntTipContactRadius )

		if ( TheContactRadius > TransitionRadius )
			//			F = 2 * E / ( 1 - v^2) * ( (indentation - ContactPoint) * TheContactRadius - sqrt(2) / pi * TheContactRadius^2 / tan(alpha) * ( pi / 2 - asin( b / TheContactRadius ) ) - TheContactRadius^3 / 3 / TipRadius + sqrt( TheContactRadius^2 - b^2) * ( sqrt(2)/ Pi * b / tan(alpha) + ( TheContactRadius^2 - b^2 ) / 3 / TipRadius ) )
			alpha *= Pi / 180			// convert to radian
			F = 2 * E / ( 1 - v^2) * ( ( Indentation - ContactPoint ) * TheContactRadius - nPyramid / pi * sin( pi / nPyramid ) * TheContactRadius^2 / 2 / tan(alpha) * ( pi / 2 - asin( TransitionRadius / TheContactRadius ) ) - TheContactRadius^3 / 3 / TipRadius + sqrt( TheContactRadius^2 - TransitionRadius^2) * ( nPyramid / Pi * sin( Pi / nPyramid) * TransitionRadius / 2 / tan(alpha) + ( TheContactRadius^2 - TransitionRadius^2 ) / 3 / TipRadius ) ) - ForceOffset
		else
			F = 4 / 3 * E / ( 1 - v^2) * sqrt(TipRadius) * (indentation - ContactPoint)^(3/2) - ForceOffset
		endif
		if (Numtype(F) == 0 )
			return(F)
		else
			return(0)
		endif
	else
		return NAN
	endif

End

//#################################################################
Function ForcevsIndHertzSphereOn3Pyramid(w,Indentation) : FitFunc
	//#################################################################
	//
	// this functions follows the paper by Felix Rico et al. Blunted Pyramidal Tips
	// PhysRev E 2005
	//
	// this paper does not give an analytic solution for F vs indent
	// instead the indentation is calculated from the contact radius
	// the force then needs the indnetation and the contact radius
	//
	// to speed up the process the indentation vs contact radius, which does not depend on E is
	// calculated before (in HertzFitDoHertzFit) and only used here as a LookUp Table
	//
	// The Rico model is usable for 3-sided and 4-sided pyramids -> nPyramid
	//
	// the transition point from sphere to pyramid is called b (distance from very end to plane of transition)
	// we assume here a smooth transition, i.e. the transition occurs were the angle of the sphere becomes identical
	// with the pyramid angle
	//
	// Manfred Radmacher Nov. 26, 2019
	//

	Wave w
	Variable Indentation

	variable E = w[0]
	variable v =  w[1]
	variable alpha = w[2]
	variable TipRadius = w[3]
	variable ContactPoint = w[4]
	variable ForceOffset = w[5]

	variable nPyramid = 3
	variable F

	// waves are initialized at the beginning of DoHertzFit
	wave BluntTipContactRadius = $(FOLDER_HERTZFIT+":BluntTipContactRadius")
	wave BluntTipIndent = $(FOLDER_HERTZFIT+":BluntTipIndent")
	NVAR TransitionRadius = $(FOLDER_HERTZFIT+":TransistionRadius")

	if (numtype(Indentation) == 0 )
		//
		//		variable TransitionRadius = CalcContactForBluntTip(TipRadius, alpha, nPyramid)

		variable TheContactRadius = interp(Indentation, BluntTipIndent, BluntTipContactRadius )
		//		print Indentation, BinarySearch(BluntTipIndent, Indentation ), BluntTipContactRadius[ BinarySearch(BluntTipIndent, Indentation )], TheContactRadius

		if ( TheContactRadius > TransitionRadius )
			//			nomenclature in Felix's paper
			//			b: contactradius at transition -> TransitionRadius
			// 			transitionradius is calculated from CalcContactForBlunt
			//			a: current contact radius -> TheContactRadius
			//				the contact raidus is read from BluntTipContactRadius looking at the same position where the current indentation is the same as BluntTipIndent

			alpha *= Pi / 180			// convert to radian
			F = 2 * E / ( 1 - v^2) * ( ( Indentation - ContactPoint ) * TheContactRadius - nPyramid / pi * sin( pi / nPyramid ) * TheContactRadius^2 / 2 / tan(alpha) * ( pi / 2 - asin( TransitionRadius / TheContactRadius ) ) - TheContactRadius^3 / 3 / TipRadius + sqrt( TheContactRadius^2 - TransitionRadius^2) * ( nPyramid / Pi * sin( Pi / nPyramid) * TransitionRadius / 2 / tan(alpha) + ( TheContactRadius^2 - TransitionRadius^2 ) / 3 / TipRadius ) ) - ForceOffset
		else
			F = 4 / 3 * E / ( 1 - v^2) * sqrt(TipRadius) * (indentation - ContactPoint)^(3/2) - ForceOffset
		endif
		if (Numtype(F) == 0 )
			return(F)
		else
			return(0)
		endif
	else
		return NAN
	endif
End

//#############################################################
Function CalcContactForBluntTip(TheTipRadius, TheTipAngle, nPyramid, MaxIndentation)
	variable TheTipRadius, TheTipAngle, nPyramid, MaxIndentation

	wave Indentation = $(FOLDER_HERTZFIT+":Indentation")
	TheTipAngle *= Pi / 180			// convert to radian
	//	variable MaxIndentation = Wavemax(Indentation)
	if ( MaxIndentation < 50e-9 )
		MaxIndentation = 50e-9
	endif
	variable MaxContactRadius = MaxIndentation * tan(TheTipAngle) * 2

	variable nPointsTestWaves = 1000
	wave BluntTipIndent = $(FOLDER_HERTZFIT+":BluntTipIndent")
	wave BluntTipContactRadius = $(FOLDER_HERTZFIT+":BluntTipContactRadius")
	redimension /n=(nPointsTestWaves) BluntTipIndent, BluntTipContactRadius

	variable DeltaContactRadius = MaxContactRadius / ( nPointsTestWaves - 1 )
	//	variable b =  TheTipRadius *  cos( TheTipAngle )		// contact radius at transition
	variable B =  TheTipRadius / tan(TheTipAngle)

	//	variable HeightAtTransition = 	b^2 / TheTipRadius		// height at contact radius, changed MR 30.3.2020
	//	variable IndentTransition = 2 * HeightAtTransition			// indentation corresponding to this situation
	variable TransitionPoint = trunc(b / DeltaContactRadius)
	if ( TransitionPoint > 	nPointsTestWaves )
		nPointsTestWaves = 2 * TransitionPoint
		make /n=(nPointsTestWaves) /O BluntTipIndent, BluntTipContactRadius
	endif
	//	if ( numtype(BluntTipContactRadius[1]) != 0 )
	BluntTipContactRadius = DeltaContactRadius * p
	BluntTipIndent = -1
	BluntTipIndent[TransitionPoint, ] = BluntTipContactRadius / tan(TheTipAngle) * nPyramid / pi * sin( pi / nPyramid ) * ( pi / 2 -  asin( b / BluntTipContactRadius ) )
	BluntTipIndent[TransitionPoint, ] -= BluntTipContactRadius / TheTipRadius * ( sqrt( BluntTipContactRadius^2 - b^2 ) - BluntTipContactRadius )

	BluntTipIndent[0, TransitionPoint] =  BluntTipContactRadius[p]^2 / TheTipRadius
	//	endif
	NVAR /Z TransitionRadius = $(FOLDER_HERTZFIT+":TransistionRadius")
	if ( NVAR_Exists(TransitionRadius) == 0 )
		Variable /G $(FOLDER_HERTZFIT+":TransistionRadius")
	endif
	TransitionRadius = b

	return(b)
end


//###########################################################################
Function HertzFitWindowModelPopUp(ctrlName,popNum,popStr) : PopupMenuControl
	//###########################################################################
	String ctrlName
	Variable popNum
	String popStr
	// 	StrConstant FitModelPopStr = "Cylinder;Sphere;Cone;3Pyramid;4Pyramid;SphereOn3Pyramid;SphereOn4Pyramid;ConeBEC;PyramidBEC;"

	NVar  WhichCurve = $(FOLDER_HERTZFIT+":WhichCurve")
	NVar  FitModel = $(FOLDER_HERTZFIT+":FitModel")
	NVAR ix = $(FOLDER_FV+":ix")
	NVAR iy = $(FOLDER_FV+":iy")
	switch (popnum)
		case 1:
			FitModel = FitModelCylinder
			break
		case 2:
			FitModel = FitModelSphere
			break
		case 3:
			FitModel = FitModelCone
			break
		case 4:
			FitModel = FitModel3Pyramid
			break
		case 5:
			FitModel = FitModel4Pyramid
			break
		case 6:
			FitModel = FitModelSphereOn3Pyramid
			break
		case 7:
			FitModel = FitModelSphereOn4Pyramid
			break
		case 8:
			FitModel = FitModelConeBEC
			break
		case 9:
			FitModel = FitModelPyramidBEC
			break
		case 10:
			FitModel = FitModelSphereBECGarcia
			break
		case 11:
			FitModel = FitModelSphereBECChadwick
			break
		default:
			FitModel = FitModel4Pyramid
			break
	endswitch
//	FVanaLoop(NO_1D_LOADING, ix, iy)
	HertzFitDoHertzfit(Approach)
	//	FVanaPrintResults("HertzFitWindowModelPopUp")
End

//###########################################################################
Function HertzFitRangeWhichDataPopUp(ctrlName,popNum,popStr) : PopupMenuControl
	//###########################################################################
	String ctrlName
	Variable popNum
	String popStr

	NVAR FitRangeWhichData = $(FOLDER_HERTZFIT+":FitRangeWhichData")
	NVAR ix = $(FOLDER_FV+":ix")
	NVAR iy = $(FOLDER_FV+":iy")

	switch (popnum)
		case 1:
			FitRangeWhichData = RangeWichDataDeflection	//deflection
			break
		case 2:
			FitRangeWhichData = RangeWichDataForce		//force
			break
		case 3:
			FitRangeWhichData = RangeWichDataIndentation //indentation
			break
		default:
			FitRangeWhichData = RangeWichDataDeflection
			break
	endswitch
	string strDeltaUnit = "m"
	if ( FitRangeWhichData == RangeWichDataForce )
		strDeltaUnit = "N"
	endif
	SetVariable setvarDelta1,format="%.1W1P"+strDeltaUnit
	SetVariable setvarDelta2,format="%.1W1P"+strDeltaUnit
	HertzFitDoHertzFit(Approach)
//	FVanaLoop(NO_1D_LOADING, ix, iy)
	//	FVanaPrintResults("HertzFitWindowFitRangeWhichDataPopUp")
End


//###########################################################################
Function HertzFitWindowWhichCurvePopUp(ctrlName,popNum,popStr) : PopupMenuControl
	//###########################################################################
	String ctrlName
	Variable popNum
	String popStr

	NVAR  WhichCurve = $(FOLDER_HERTZFIT+":WhichCurve")
//	NVAR ix = $(FOLDER_FV+":ix")
//	NVAR iy = $(FOLDER_FV+":iy")
	WhichCurve = popNum
	HertzFitDoHertzFit(Approach)
//	FVanaLoop(NO_1D_LOADING, ix, iy)
	//	FVanaPrintResults("HertzFitWindowWhichCurvePopUp")
End


//===============================================
Function CheckProcIterateFit(ctrlName,checked) : CheckBoxControl
	//===============================================
	String ctrlName
	Variable checked
//	NVAR ix = $(FOLDER_FV+":ix")
//	NVAR iy = $(FOLDER_FV+":iy")
	HertzFitDoHertzFit(Approach)
//	FVanaLoop(NO_1D_LOADING, ix, iy)
	//	FVanaPrintResults("CheckProcIterateFit")
End


//===============================================
Function CheckProcTiltCorr(ctrlName,checked) : CheckBoxControl
	//===============================================
	String ctrlName
	Variable checked
//	NVAR ix = $(FOLDER_FV+":ix")
//	NVAR iy = $(FOLDER_FV+":iy")
	HertzFitDoHertzFit(Approach)
//	FVanaLoop(NO_1D_LOADING, ix, iy)
	//	FVanaPrintResults("CheckProcTiltCorr")
End


//===============================================
Function CheckProcPrintResults(ctrlName,checked) : CheckBoxControl
	//===============================================
	String ctrlName
	Variable checked
//	NVAR ix = $(FOLDER_FV+":ix")
//	NVAR iy = $(FOLDER_FV+":iy")
	HertzFitDoHertzFit(Approach)
//	FVanaLoop(NO_1D_LOADING, ix, iy)
	//	FVanaPrintResults("CheckProcPrintResults")
End
