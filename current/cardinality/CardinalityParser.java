package cardinality;

import parse.Parse;
import parse.Parser;

/**
 * Centralized parser for cardinality estimator flags.
 * <p>
 * Consolidates flag parsing duplicated across DDLCalibrationDriver2 and
 * LowComplexityCalibrationDriver into a single location.
 *
 * @author Nahida, Brian Bushnell
 * @date April 2026
 */
public class CardinalityParser {

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parse a cardinality-related flag.
	 * @param arg raw argument string (unused, retained for caller compatibility)
	 * @param a flag name (lowercase)
	 * @param b flag value string
	 * @return true if the flag was recognized and handled
	 */
	public static boolean parse(final String arg, final String a, final String b){
		if(a.equals("loglogtype") || a.equals("type")){Parser.loglogType=b.toLowerCase();}
		else if(a.equals("cf") || a.equals("loglogcf")){CorrectionFactor.USE_CORRECTION=Parse.parseBoolean(b);}
		else if(a.equals("cardcf")){CardinalityTracker.USE_CARD_CF=Parse.parseBoolean(b);}
		else if(a.equals("dlccf")){CorrectionFactor.USE_DLC_CF=Parse.parseBoolean(b);}
		else if(a.equals("printcv") || a.equals("cv")){DDLCalibrationDriver.PRINT_CV=Parse.parseBoolean(b);}
		else if(a.equals("printstd") || a.equals("std")){DDLCalibrationDriver.PRINT_STD=Parse.parseBoolean(b);}
		else if(a.equals("formulas") || a.equals("useformulas")){CorrectionFactor.USE_FORMULAS=Parse.parseBoolean(b);}
		else if(a.equals("sbsformula") || a.equals("usesbsformula")){CorrectionFactor.USE_SBS_FORMULA=Parse.parseBoolean(b);}
		else if(a.equals("meancfformula") || a.equals("usemeancfformula")){CorrectionFactor.USE_MEAN_CF_FORMULA=Parse.parseBoolean(b);}
		else if(a.equals("hccfformula") || a.equals("usehccfformula")){CorrectionFactor.USE_HC_CF_FORMULA=Parse.parseBoolean(b);}
		else if(a.equals("tracecf")){CorrectionFactor.TRACE_CF=Parse.parseBoolean(b);}
		else if(a.equals("lchistfile")){CorrectionFactor.sbsFile=b;}
		else if(a.equals("lchistmultfile")){CorrectionFactor.sbsMultFile=b;}
		else if(a.equals("histlc") || a.equals("historylc")){CardinalityTracker.USE_HISTORY_FOR_LC=Parse.parseBoolean(b);}
		else if(a.equals("lchisthybrid") || a.equals("lchybrid") || a.equals("usesbs") || a.equals("sbsinhybrid")){CardinalityTracker.USE_SBS_IN_HYBRID=Parse.parseBoolean(b);}
		else if(a.equals("hllhistcf") || a.equals("histcf")){CardinalityStats.HLL_HIST_TERMINAL_CF=Double.parseDouble(b);}
		else if(a.equals("tmcf") || a.equals("terminalmeancf")){AbstractCardStats.OVERRIDE_TERMINAL_MEAN_CF=Float.parseFloat(b);}
		else if(a.equals("tmpcf") || a.equals("terminalmeanpluscf")){AbstractCardStats.OVERRIDE_TERMINAL_MEAN_PLUS_CF=Float.parseFloat(b);}
		else if(a.equals("microlc")){CardinalityStats.USE_MICRO_FOR_LC=Parse.parseBoolean(b); AbstractCardStats.USE_MICRO_FOR_LC=Parse.parseBoolean(b);}
		else if(a.equals("usemicro") || a.equals("usemicroindex")){AbstractCardStats.USE_MICRO_INDEX=Parse.parseBoolean(b);}

		// DLC/CardStats tuning flags
		else if(a.equals("dlcalpha") || a.equals("alpha")){CardinalityStats.DLC_ALPHA=Float.parseFloat(b); AbstractCardStats.DLC_ALPHA=Float.parseFloat(b);}
		else if(a.equals("dlctarget")){CardinalityStats.DLC_TARGET_FRAC=Float.parseFloat(b); AbstractCardStats.DLC_TARGET_FRAC=Float.parseFloat(b);}
		else if(a.equals("dlcblendlo")){CardinalityStats.DLC_BLEND_LO=Float.parseFloat(b); AbstractCardStats.DLC_BLEND_LO=Float.parseFloat(b);}
		else if(a.equals("dlcblendhi")){CardinalityStats.DLC_BLEND_HI=Float.parseFloat(b); AbstractCardStats.DLC_BLEND_HI=Float.parseFloat(b);}
		else if(a.equals("cfiters") || a.equals("cfiterations")){CardinalityStats.DEFAULT_CF_ITERS=Integer.parseInt(b); AbstractCardStats.DEFAULT_CF_ITERS=Integer.parseInt(b);}
		else if(a.equals("cfdif") || a.equals("cfconvergence")){CardinalityStats.DEFAULT_CF_DIF=Double.parseDouble(b); AbstractCardStats.DEFAULT_CF_DIF=Double.parseDouble(b);}
		else if(a.equals("cfmult") || a.equals("minseedmult")){CardinalityStats.MIN_SEED_CF_MULT=Float.parseFloat(b); AbstractCardStats.MIN_SEED_CF_MULT=Float.parseFloat(b);}
		else if(a.equals("minvfraction") || a.equals("minvk")){
			final float x=Float.parseFloat(b);
			if(x<1){CardinalityStats.DLC_MIN_VK_FRACTION=x; AbstractCardStats.DLC_MIN_VK_FRACTION=x;}
			else{CardinalityStats.DLC_MIN_VK=(int)x; AbstractCardStats.DLC_MIN_VK=(int)x;}
		}
		else if(a.equals("dlcinfopow") || a.equals("infopow")){AbstractCardStats.DLC_INFO_POWER=Float.parseFloat(b);}
		else if(a.equals("hcinfopow") || a.equals("hcpow")){AbstractCardStats.HC_INFO_POWER=Float.parseFloat(b);}
		else if(a.equals("dlchistblendlo")){AbstractCardStats.DLCSBS_BLEND_LO=Float.parseFloat(b);}
		else if(a.equals("dlchistblendhi")){AbstractCardStats.DLCSBS_BLEND_HI=Float.parseFloat(b);}
		else if(a.equals("hybridblendlo") || a.equals("hblo")){AbstractCardStats.HYBRID_BLEND_LO=Float.parseFloat(b);}
		else if(a.equals("hybridblendhi") || a.equals("hbhi")){AbstractCardStats.HYBRID_BLEND_HI=Float.parseFloat(b);}
		else if(a.equals("hybridblendlog") || a.equals("hblog")){AbstractCardStats.HYBRID_BLEND_LOG=Parse.parseBoolean(b);}
		else if(a.equals("ldlcblo")){AbstractCardStats.LDLC_B_LO=Float.parseFloat(b);}
		else if(a.equals("ldlcbhi")){AbstractCardStats.LDLC_B_HI=Float.parseFloat(b);}
		else if(a.equals("dlcinfomode") || a.equals("infomode")){AbstractCardStats.DLC_INFO_MODE=Integer.parseInt(b);}

		// Per-class flags
		else if(a.equals("earlypromote") || a.equals("ep")){
			DynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
			DynamicLogLog4.EARLY_PROMOTE=Parse.parseBoolean(b);
			BankedDynamicLogLog3.EARLY_PROMOTE=Parse.parseBoolean(b);
			DynamicLogLog3v4.EARLY_PROMOTE=Parse.parseBoolean(b);
		}else if(a.equals("promotethreshold") || a.equals("pt")){
			DynamicLogLog3.PROMOTE_THRESHOLD=Integer.parseInt(b);
			DynamicLogLog4.PROMOTE_THRESHOLD=Integer.parseInt(b);
			DynamicLogLog3v4.PROMOTE_THRESHOLD=Integer.parseInt(b);
		}else if(a.equals("promotefrac") || a.equals("pf")){DynamicLogLog3v4.PROMOTE_FRAC=BankedDynamicLogLog3.PROMOTE_FRAC=BankedDynamicLogLog4.PROMOTE_FRAC=BankedDynamicLogLog5.PROMOTE_FRAC=Float.parseFloat(b);}
		else if(a.equals("overflowscale") || a.equals("os")){
			DynamicLogLog3.OVERFLOW_SCALE=Double.parseDouble(b);
			BankedDynamicLogLog3.OVERFLOW_SCALE=Double.parseDouble(b);
			DynamicLogLog3v4.OVERFLOW_SCALE=Double.parseDouble(b);
		}else if(a.equals("usestoredoverflow") || a.equals("uso")){
			DynamicLogLog3.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
			BankedDynamicLogLog3.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
			DynamicLogLog3v4.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
		}
		else if(a.equals("edllhbits") || a.equals("ehbits")){ExpandedDynamicLogLog8.EDLL_HBITS=Integer.parseInt(b); ExpandedDynamicLogLog9.EDLL_HBITS=Integer.parseInt(b);}
		else if(a.equals("etllhsb") || a.equals("hsbmode")){StateTable.ETLL_HSB_MODE=Integer.parseInt(b);}
		else if(a.equals("hsbetll3")){
			final String[] parts=b.split(",");
			final double[] tbl=new double[parts.length];
			for(int i=0; i<parts.length; i++){tbl[i]=Double.parseDouble(parts[i]);}
			StateTable.CF_ETLL_3_OVERRIDE=tbl;
		}
		else if(a.equals("calcfgra") || a.equals("fgra")){UltraDynamicLogLog6.CALC_FGRA=Parse.parseBoolean(b);}
		else if(a.equals("saturate") || a.equals("sat")){UltraDynamicLogLog6.SATURATE_ON_OVERFLOW=Parse.parseBoolean(b);}
		else if(a.equals("printdlctiers")){DDLCalibrationDriver.PRINT_DLC_TIERS=Parse.parseBoolean(b);}
		else if(a.equals("capdenom") || a.equals("cd")){HalfCappedDynamicLogLog4.CAP_DENOM=Integer.parseInt(b);}
		else if(a.equals("demotionmode") || a.equals("dm")){VariableCompressedDynamicLogLog4.DEMOTION_MODE=Integer.parseInt(b); ArithmeticVariableDynamicLogLog32.DEMOTION_MODE=Integer.parseInt(b); ArithmeticVariableDynamicLogLog36.DEMOTION_MODE=Integer.parseInt(b); ArithmeticVariableDynamicLogLog34.DEMOTION_MODE=Integer.parseInt(b);}
		// else if(a.equals("vcdll4histlimit") || a.equals("vhl")){VariableCompressedDynamicLogLog4.HIST_LIMIT=Integer.parseInt(b);}
		else if(a.equals("histtiers") || a.equals("htl")){ArithmeticVariableDynamicLogLog32.HIST_TIER_LIMIT=Integer.parseInt(b); ArithmeticVariableDynamicLogLog32.reconfigure(); ArithmeticVariableDynamicLogLog36.HIST_TIER_LIMIT=Integer.parseInt(b); ArithmeticVariableDynamicLogLog36.reconfigure(); ArithmeticVariableDynamicLogLog34.HIST_TIER_LIMIT=Integer.parseInt(b); ArithmeticVariableDynamicLogLog34.reconfigure();}
		else if(a.equals("numtails")){CompressedTwinTailLogLog.NUM_TAILS=Integer.parseInt(b); CompressedTwinTailLogLog.reconfigure();}
		else if(a.equals("histlen")){CompressedTwinTailLogLog.HIST_LEN=Integer.parseInt(b); CompressedTwinTailLogLog.reconfigure();}
		else if(a.equals("expanded")){
			CompressedTwinTailLogLog.TIER_NUMER=2; CompressedTwinTailLogLog.TIER_DENOM=1;
			CompressedTwinTailLogLog.EXP_BITS=5; CompressedTwinTailLogLog.reconfigure();
		}

		// PLL family flags
		else if(a.equals("plloffset") || a.equals("pco")){
			ProtoLogLog16c.CF_OFFSET=Double.parseDouble(b);
		}else if(a.equals("hbits")){final int v=Integer.parseInt(b); ProtoLogLog16c.HISTORY_BITS=v;}
		else if(a.equals("lbits")){final int v=Integer.parseInt(b); ProtoLogLog16c.LUCK_BITS=v;}
		else if(a.equals("mbits")){final int v=Integer.parseInt(b); ProtoLogLog16c.MANTISSA_BITS=v;}
		else if(a.equals("abits")){final int v=Integer.parseInt(b); ProtoLogLog16c.ANDTISSA_BITS=v;}
		else if(a.equals("nbits")){final int v=Integer.parseInt(b); ProtoLogLog16c.NLZ2_BITS=v;}

		// DynamicDemiLog8 flags
		else if(a.equals("empiricalmantissa") || a.equals("em")){DynamicDemiLog8.USE_EMPIRICAL_MANTISSA=Parse.parseBoolean(b);}
		else if(a.equals("mantissacfoffset") || a.equals("mco")){DynamicDemiLog8.MANTISSA_CF_OFFSET=Double.parseDouble(b);}

		// FLL2 flags
		else if(a.equals("fll2mult")){FutureLogLog2.TERMINAL_CORRECTION=Double.parseDouble(b);}
		else if(a.equals("clampoverflow")){FutureLogLog2.CLAMP_OVERFLOW=Parse.parseBoolean(b);}
		else if(a.equals("hcweight") || a.equals("ldlcweight")){CardinalityTracker.OVERRIDE_LDLC_HC_WEIGHT=Double.parseDouble(b);}
		else if(a.equals("hybridweight") || a.equals("hw")){CardinalityTracker.OVERRIDE_HLDLC_WEIGHT=Float.parseFloat(b);}
		else if(a.equals("hcscale")){AbstractCardStats.HC_SCALE=Float.parseFloat(b); AbstractCardStats.HC_SCALE_EXPLICIT=true;}

		// BCLL flags
		else if(a.equals("bcllmult")){BankedCeilingLogLog.TERMINAL_CORRECTION=Double.parseDouble(b);}
		else if(a.equals("cvpow") || a.equals("cvpower")){BankedCeilingLogLog.CV_POWER=Double.parseDouble(b); DynamicLogLog4.CV_POWER=Double.parseDouble(b); TwinTailLogLog.CV_POWER=Double.parseDouble(b);}

		// TTLL flags
		else if(a.equals("ttllmult")){TwinTailLogLog.TERMINAL_CORRECTION=Double.parseDouble(b);}
		else if(a.equals("ttllcvpow")){TwinTailLogLog.CV_POWER=Double.parseDouble(b);}

		// DLL4 word table flags
		else if(a.equals("dll4wordmult")){DynamicLogLog4.WORD_TERMINAL_CORRECTION=Double.parseDouble(b);}
		else if(a.equals("dll4cvpow")){DynamicLogLog4.CV_POWER=Double.parseDouble(b);}

		// Cascading overflow flags (co/io are mutually exclusive; "co" alias is handled by correctoverflow)
		else if(a.equals("correctoverflow") || a.equals("co")){
			final boolean x=Parse.parseBoolean(b);
			DynamicLogLog3.CORRECT_OVERFLOW=x;
			BankedDynamicLogLog3.CORRECT_OVERFLOW=x;
			BankedDynamicLogLog4.CORRECT_OVERFLOW=x;
			CompressedDynamicLogLog3.CORRECT_OVERFLOW=x;
			DynamicLogLog3v4.CORRECT_OVERFLOW=x;
			if(x){DynamicLogLog3.IGNORE_OVERFLOW=DynamicLogLog2.IGNORE_OVERFLOW=BankedDynamicLogLog3.IGNORE_OVERFLOW=BankedDynamicLogLog4.IGNORE_OVERFLOW=DynamicLogLog3v4.IGNORE_OVERFLOW=false;}
		}else if(a.equals("ignoreoverflow") || a.equals("io")){
			final boolean x=Parse.parseBoolean(b);
			DynamicLogLog3.IGNORE_OVERFLOW=x;
			DynamicLogLog2.IGNORE_OVERFLOW=x;
			BankedDynamicLogLog3.IGNORE_OVERFLOW=x;
			BankedDynamicLogLog4.IGNORE_OVERFLOW=x;
			CompressedDynamicLogLog3.IGNORE_OVERFLOW=x;
			DynamicLogLog3v4.IGNORE_OVERFLOW=x;
			if(x){DynamicLogLog3.CORRECT_OVERFLOW=BankedDynamicLogLog3.CORRECT_OVERFLOW=BankedDynamicLogLog4.CORRECT_OVERFLOW=DynamicLogLog3v4.CORRECT_OVERFLOW=false;}
		}else if(a.equals("iobias") || a.equals("iob")){
			final boolean x=Parse.parseBoolean(b);
			DynamicLogLog3.USE_IO_BIAS=x;
			DynamicLogLog3v4.USE_IO_BIAS=x;
			BankedDynamicLogLog3.USE_IO_BIAS=x;
		}else if(a.equals("iobfile") || a.equals("iobf")){
			DynamicLogLog3.loadIOBias(b);
			DynamicLogLog3v4.loadIOBias(b);
			BankedDynamicLogLog3.loadIOBias(b);
		}else if(a.equals("clamptoadded") || a.equals("clamp")){CardinalityTracker.clampToAdded=Parse.parseBoolean(b);}
		else if(a.equals("dual")){
//			CompressedDynamicLogLog3.DUAL=Parse.parseBoolean(b);//Do nothing
		}else{return false;}
		return true;
	}

	/**
	 * Initialize global cardinality state after parsing is complete.
	 * <p>
	 * Sets up PLL mode, HC weight defaults, per-class formula coefficients,
	 * loads CF/SBS tables, and publishes the immutable CF snapshot for worker threads.
	 *
	 * @param loglogtype estimator type string (e.g. "dll4", "dll3", "bdll3")
	 * @param buckets number of buckets
	 * @param k hash k-mer length
	 * @param cffile correction factor file path, or null for default
	 * @param pllmode PLL mode string, or null to skip
	 */
	public static void initializeAll(final String loglogtype, final int buckets, final int k,
			String cffile, final String pllmode){

		// PLL mode application
		if(pllmode!=null){
			if(pllmode.equals("mantissa")){ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_MANTISSA);}
			else if(pllmode.equals("andtissa")){ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_ANDTISSA);}
			else if(pllmode.equals("nlz2")){ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_NLZ2);}
			else if(pllmode.equals("history")){ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_HISTORY);}
			else if(pllmode.equals("luck")){ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_LUCK);}
			else if(pllmode.equals("none")){ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_NONE);}
			else if(pllmode.equals("histmant") || pllmode.equals("historymantissa")){ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_HISTORY|ProtoLogLog16c.MODE_MANTISSA);}
			else{throw new RuntimeException("Unknown pllmode: "+pllmode);}
		}

		// Set LDLC_HC_WEIGHT from the tracker's instance method (or command-line override)
		{
			final CardinalityTracker probe=CardinalityTracker.makeTracker(loglogtype, 64, 31, 0, 0);
			CardinalityTracker.LDLC_HC_WEIGHT=(CardinalityTracker.OVERRIDE_LDLC_HC_WEIGHT>=0
				? CardinalityTracker.OVERRIDE_LDLC_HC_WEIGHT
				: probe.ldlcHcWeight());
		}

		// Class init and column whitelist
		DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, 0L, 0);
		DDLCalibrationDriver.v3ColsForType(loglogtype);

		// Reset per-class corrections (preserve explicit command-line overrides)
		if(!AbstractCardStats.HC_SCALE_EXPLICIT){AbstractCardStats.HC_SCALE=1.0f;}
		StateTable.USE_AUDLL32_HSB=false;

		// Per-class Mean CF formula coefficients
		boolean classHasMeanCf=true, classHasHcCf=false, classHasHmeanmCf=false;
		if(loglogtype.equals("cdll4")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.5; CorrectionFactor.sbsFile=CompressedDynamicLogLog4.SBS_FILE;}
		else if(loglogtype.equals("cdll5")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.5;
			CorrectionFactor.sbsFile=CompressedDynamicLogLog5.SBS_FILE;
		}else if(loglogtype.equals("bcdll5")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.5;
			CorrectionFactor.sbsFile=BankedCompressedDynamicLogLog5.SBS_FILE;
		}else if(loglogtype.equals("acdll5")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.5;
			CorrectionFactor.sbsFile=ArithmeticCompressedDynamicLogLog5.SBS_FILE;
		}else if(loglogtype.equals("vcdll4")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.5; CorrectionFactor.sbsFile=VariableCompressedDynamicLogLog4.SBS_FILE;}
		else if(loglogtype.equals("acdll4")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.5;
			CorrectionFactor.sbsFile=ArithmeticCompressedDynamicLogLog4.SBS_FILE;
		}else if(loglogtype.equals("hcdll4")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.5;
			if(HalfCappedDynamicLogLog4.SBS_FILE!=null){CorrectionFactor.sbsFile=HalfCappedDynamicLogLog4.SBS_FILE;}
		}else if(loglogtype.equals("edll8")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=0.5;
			if(CorrectionFactor.sbsFile==null && ExpandedDynamicLogLog8.SBS_FILE!=null){CorrectionFactor.sbsFile=ExpandedDynamicLogLog8.SBS_FILE;}
		}else if(loglogtype.equals("edll9")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=0.5;
			if(CorrectionFactor.sbsFile==null && ExpandedDynamicLogLog9.SBS_FILE!=null){CorrectionFactor.sbsFile=ExpandedDynamicLogLog9.SBS_FILE;}
			if(!AbstractCardStats.HC_SCALE_EXPLICIT){AbstractCardStats.HC_SCALE=1.030880f;}
		}else if(loglogtype.equals("audll32")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.0;
			StateTable.USE_AUDLL32_HSB=true;
			CorrectionFactor.sbsFile=ArithmeticUltraDynamicLogLog32.SBS_FILE;
			if(!AbstractCardStats.HC_SCALE_EXPLICIT){AbstractCardStats.HC_SCALE=1.008f;}
		}else if(loglogtype.equals("audll33")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.0;
			StateTable.USE_AUDLL32_HSB=true;
			CorrectionFactor.sbsFile=ArithmeticUltraDynamicLogLog32.SBS_FILE;
			if(!AbstractCardStats.HC_SCALE_EXPLICIT){AbstractCardStats.HC_SCALE=1.008f;}
		}else if(loglogtype.equals("avdll32")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_UDLL6;
			CorrectionFactor.meanhCfCoeffs=CorrectionFactor.MCF_UDLL6_MEANH;
			classHasHcCf=true;
		}else if(loglogtype.equals("avdll36")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_UDLL6;
			CorrectionFactor.meanhCfCoeffs=CorrectionFactor.MCF_UDLL6_MEANH;
			classHasHcCf=true;
		}else if(loglogtype.equals("avdll34")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_UDLL6;
			CorrectionFactor.meanhCfCoeffs=CorrectionFactor.MCF_UDLL6_MEANH;
			classHasHcCf=true;
		}else if(loglogtype.equals("dhdll3") || loglogtype.equals("cdll3")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=(CompressedDynamicLogLog3.DUAL ? 2 : 1.5);}
		else if(loglogtype.equals("bcdll3")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.5;}
		else if(loglogtype.equals("cttll")){AbstractCardStats.TIER_SCALE=CompressedTwinTailLogLog.TIER_SCALE; CorrectionFactor.sbsFile=CompressedTwinTailLogLog.sbsFile();}
		else if(loglogtype.equals("ettll5")){AbstractCardStats.TIER_SCALE=ExpandedTwinTailLogLog5.TIER_SCALE; CorrectionFactor.sbsFile=ExpandedTwinTailLogLog5.SBS_FILE;}
		else if(loglogtype.equals("dhdll4")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=2;}
		else if(loglogtype.equals("dll4") || loglogtype.equals("dll4m")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4;}
		else if(loglogtype.equals("ll6")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_LL6;}
		else if(loglogtype.equals("dll3")){
			if(DynamicLogLog3.IGNORE_OVERFLOW){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL3_IOF;}
			else{CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL3_IOT;}
		}else if(loglogtype.equals("dll3v4")){
			if(DynamicLogLog3v4.IGNORE_OVERFLOW){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL3_IOF;}
			else{CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL3_IOT;}
		}else if(loglogtype.equals("dll3v2") || loglogtype.equals("dll3v3")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL3_IOT;
		}else if(loglogtype.equals("dll2")){
			if(DynamicLogLog2.IGNORE_OVERFLOW){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL2_IOF;}
			else{CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL2_IOT;}
		}else if(loglogtype.equals("bdll3")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_BDLL3_COF;}
		else if(loglogtype.equals("bdll4")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_BDLL3_COF;}
		else if(loglogtype.equals("bdll5")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_BDLL3_COF; CorrectionFactor.sbsFile=BankedDynamicLogLog5.SBS_FILE;}
		else if(loglogtype.equals("udll6") || loglogtype.equals("udll6m") || loglogtype.equals("udll36")){
			CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_UDLL6;
			CorrectionFactor.meanhCfCoeffs=CorrectionFactor.MCF_UDLL6_MEANH;
			classHasHcCf=true;
		}else if(loglogtype.equals("ddl") || loglogtype.equals("ddl10")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DDL; CorrectionFactor.hmeanmCfCoeffs=CorrectionFactor.HMCF_DDL; classHasHmeanmCf=true;}
		else if(loglogtype.equals("ddl8")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DDL8; CorrectionFactor.hmeanmCfCoeffs=CorrectionFactor.HMCF_DDL8; classHasHmeanmCf=true;}
		else{classHasMeanCf=false;}
		if(CorrectionFactor.USE_FORMULAS){
			if(!classHasMeanCf){CorrectionFactor.USE_MEAN_CF_FORMULA=false;}
			if(!classHasHcCf){CorrectionFactor.USE_HC_CF_FORMULA=false;}
			if(!classHasHmeanmCf){CorrectionFactor.USE_HMEANM_CF_FORMULA=false;}
		}

		// Load SBS, FLL2, BCLL, and DLL4 word tables
		CorrectionFactor.loadSbsTable();
		CorrectionFactor.loadSbsMultTable();
		if("fll2".equals(loglogtype)){FutureLogLog2.loadCFTable(); FutureLogLog2.loadCardCFTable();}
		if("bcll".equals(loglogtype)){BankedCeilingLogLog.loadCFTable(); BankedCeilingLogLog.loadCardCFTable();}
		if("ttll".equals(loglogtype)){TwinTailLogLog.loadCFTable(); TwinTailLogLog.loadCardCFTable();}
		if("cttll".equals(loglogtype)){CompressedTwinTailLogLog.loadPerTierTable(); CompressedTwinTailLogLog.computeHSB();}
		if("ettll5".equals(loglogtype)){ExpandedTwinTailLogLog5.loadPerTierTable();}
		if("dll4".equals(loglogtype) || "dll4m".equals(loglogtype)){DynamicLogLog4.loadWordTable();}

		// Load CF table: explicit cffile overrides, otherwise auto-select per type+mode
		if(cffile==null){cffile=defaultCFFile(loglogtype);}
		if(cffile!=null){
			CorrectionFactor.initialize(cffile, buckets);
			if("pll16c".equals(loglogtype)){ProtoLogLog16c.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("udll6".equals(loglogtype)){UltraDynamicLogLog6.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
				else if("udll6m".equals(loglogtype)){UltraDynamicLogLog6m.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
				else if("udll36".equals(loglogtype)){UltraDynamicLogLog36.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("bdll4".equals(loglogtype)){BankedDynamicLogLog4.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("bdll5".equals(loglogtype)){BankedDynamicLogLog5.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("cdll5".equals(loglogtype)){CompressedDynamicLogLog5.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("bcdll5".equals(loglogtype)){BankedCompressedDynamicLogLog5.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("acdll5".equals(loglogtype)){ArithmeticCompressedDynamicLogLog5.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("vcdll4".equals(loglogtype)){VariableCompressedDynamicLogLog4.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("acdll4".equals(loglogtype)){ArithmeticCompressedDynamicLogLog4.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("hcdll4".equals(loglogtype)){HalfCappedDynamicLogLog4.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("edll8".equals(loglogtype)){ExpandedDynamicLogLog8.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("edll9".equals(loglogtype)){ExpandedDynamicLogLog9.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("audll32".equals(loglogtype)){ArithmeticUltraDynamicLogLog32.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("audll33".equals(loglogtype)){ArithmeticUltraDynamicLogLog33.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("avdll32".equals(loglogtype)){ArithmeticVariableDynamicLogLog32.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("avdll36".equals(loglogtype)){ArithmeticVariableDynamicLogLog36.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("avdll34".equals(loglogtype)){ArithmeticVariableDynamicLogLog34.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("ertl".equals(loglogtype)){ErtlULL.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
		}

		// Publish immutable CF snapshot for worker threads
		CorrectionFactor.publishSnapshot();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Select the default CF table file for a given estimator type and overflow mode.
	 * Returns null if no table exists for the type.
	 */
	private static String defaultCFFile(final String loglogtype){
		if(loglogtype.equals("cdll4")){
			return CompressedDynamicLogLog4.CF_FILE;
		}else if(loglogtype.equals("cdll5")){
			return CompressedDynamicLogLog5.CF_FILE;
		}else if(loglogtype.equals("bcdll5")){
			return BankedCompressedDynamicLogLog5.CF_FILE;
		}else if(loglogtype.equals("acdll5")){
			return ArithmeticCompressedDynamicLogLog5.CF_FILE;
		}else if(loglogtype.equals("acdll4")){
			return ArithmeticCompressedDynamicLogLog4.CF_FILE;
		}else if(loglogtype.equals("hcdll4")){
			return HalfCappedDynamicLogLog4.CF_FILE;
		}else if(loglogtype.equals("vcdll4")){
			return VariableCompressedDynamicLogLog4.CF_FILE;
		// }else if(loglogtype.equals("edll8")){
		//	return ExpandedDynamicLogLog8.CF_FILE;
		}else if(loglogtype.equals("edll9")){
			return ExpandedDynamicLogLog9.CF_FILE;
		}else if(loglogtype.equals("audll32")){
			return ArithmeticUltraDynamicLogLog32.CF_FILE;
		}else if(loglogtype.equals("audll33")){
			return ArithmeticUltraDynamicLogLog33.CF_FILE;
		}else if(loglogtype.equals("avdll32")){
			return ArithmeticVariableDynamicLogLog32.CF_FILE;
		}else if(loglogtype.equals("avdll36")){
			return ArithmeticVariableDynamicLogLog36.CF_FILE;
		}else if(loglogtype.equals("avdll34")){
			return ArithmeticVariableDynamicLogLog34.CF_FILE;
		}else if(loglogtype.equals("dhdll3") || loglogtype.equals("cdll3") || loglogtype.equals("dhdll4")){
			return CompressedDynamicLogLog3.CF_FILE;
		}else if(loglogtype.equals("bcdll3")){
			return BankedCompressedDynamicLogLog3.CF_FILE;
		}else if(loglogtype.equals("dll3")){
			if(DynamicLogLog3.IGNORE_OVERFLOW){return "?cardinalityCorrectionDLL3_iot.tsv.gz";}
			else if(!DynamicLogLog3.CORRECT_OVERFLOW){return "?cardinalityCorrectionDLL3_cof.tsv.gz";}
			else{return "?cardinalityCorrectionDLL3_cot.tsv.gz";}
		}else if(loglogtype.equals("dll4") || loglogtype.equals("dll4m")){
			return DynamicLogLog4.CF_FILE;
		}else if(loglogtype.equals("dll2")){
			if(DynamicLogLog2.IGNORE_OVERFLOW){return "?cardinalityCorrectionDLL2_iot.tsv.gz";}
			else{return "?cardinalityCorrectionDLL2_iof.tsv.gz";}
		}else if(loglogtype.equals("dll3v4")){
			if(DynamicLogLog3v4.IGNORE_OVERFLOW){return "?cardinalityCorrectionDLL3v4_iot.tsv.gz";}
			else if(!DynamicLogLog3v4.CORRECT_OVERFLOW){return "?cardinalityCorrectionDLL3v4_cof.tsv.gz";}
			else{return "?cardinalityCorrectionDLL3v4_cot.tsv.gz";}
		}else if(loglogtype.equals("dll3v2") || loglogtype.equals("dll3v3")){
			return DynamicLogLog3.CF_FILE;
		}else if(loglogtype.equals("bdll3")){
			if(BankedDynamicLogLog3.IGNORE_OVERFLOW){return "?cardinalityCorrectionBDLL3_iot.tsv.gz";}
			else if(!BankedDynamicLogLog3.CORRECT_OVERFLOW){return "?cardinalityCorrectionBDLL3_cof.tsv.gz";}
			else{return "?cardinalityCorrectionBDLL3_cot.tsv.gz";}
		}else if(loglogtype.equals("bdll4")){
			return BankedDynamicLogLog4.CF_FILE;
		}else if(loglogtype.equals("bdll5")){
			return BankedDynamicLogLog5.CF_FILE;
		}else if(loglogtype.equals("ll6")){
			return LogLog6.CF_FILE;
		}else if(loglogtype.equals("udll6") || loglogtype.equals("udll6m") || loglogtype.equals("udll36")){
			return UltraDynamicLogLog6.CF_FILE;
		}else if(loglogtype.equals("ddl") || loglogtype.equals("ddl10")){
			return DynamicDemiLog.CF_FILE;
		}else if(loglogtype.equals("ddl8")){
			return DynamicDemiLog8.CF_FILE;
		}else if(loglogtype.equals("ttll")){
			return "?"+TwinTailLogLog.CF_FILE;
		}else if(loglogtype.equals("ttll3")){
			return "?cardinalityCorrectionTTLL3.tsv.gz";
		}else if(loglogtype.equals("ttll4")){
			return "?cardinalityCorrectionTTLL4.tsv.gz";
		}else if(loglogtype.equals("cttll")){
			if(CompressedTwinTailLogLog.NUM_TAILS==4){return "?cardinalityCorrectionCTTLL_quad.tsv.gz";}
			return "?cardinalityCorrectionCTTLL.tsv.gz";
		}else if(loglogtype.equals("ettll5")){
			if(!CorrectionFactor.USE_CORRECTION){return null;}
			return "?cardinalityCorrectionETTLL5.tsv.gz";
		}
		return null;
	}

}
