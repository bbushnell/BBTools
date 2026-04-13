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

	/**
	 * Parse a cardinality-related flag.
	 * @param a flag name (lowercase)
	 * @param b flag value string
	 * @return true if the flag was recognized and handled
	 */
	public static boolean parse(String arg, String a, String b){
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
			float x=Float.parseFloat(b);
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
			DynamicLogLog3v3.EARLY_PROMOTE=Parse.parseBoolean(b);
			DynamicLogLog3v4.EARLY_PROMOTE=Parse.parseBoolean(b);
		}
		else if(a.equals("promotethreshold") || a.equals("pt")){
			DynamicLogLog3.PROMOTE_THRESHOLD=Integer.parseInt(b);
			DynamicLogLog3v2.PROMOTE_THRESHOLD=Integer.parseInt(b);
			DynamicLogLog4.PROMOTE_THRESHOLD=Integer.parseInt(b);
			DynamicLogLog3v3.PROMOTE_THRESHOLD=Integer.parseInt(b);
			DynamicLogLog3v4.PROMOTE_THRESHOLD=Integer.parseInt(b);
		}
		else if(a.equals("promotefrac") || a.equals("pf")){DynamicLogLog3v2.PROMOTE_FRAC=Float.parseFloat(b); DynamicLogLog3v3.PROMOTE_FRAC=Float.parseFloat(b); DynamicLogLog3v4.PROMOTE_FRAC=Float.parseFloat(b); BankedDynamicLogLog3.PROMOTE_FRAC=Float.parseFloat(b);}
		else if(a.equals("resetonpromote") || a.equals("rop")){DynamicLogLog3v2.RESET_ON_PROMOTE=Parse.parseBoolean(b);}
		else if(a.equals("overflowscale") || a.equals("os")){
			DynamicLogLog3.OVERFLOW_SCALE=Double.parseDouble(b);
			BankedDynamicLogLog3.OVERFLOW_SCALE=Double.parseDouble(b);
			DynamicLogLog3v3.OVERFLOW_SCALE=Double.parseDouble(b);
			DynamicLogLog3v4.OVERFLOW_SCALE=Double.parseDouble(b);
		}
		else if(a.equals("usestoredoverflow") || a.equals("uso")){
			DynamicLogLog3.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
			BankedDynamicLogLog3.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
			DynamicLogLog3v3.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
			DynamicLogLog3v4.USE_STORED_OVERFLOW=Parse.parseBoolean(b);
		}
		else if(a.equals("calcfgra") || a.equals("fgra")){UltraDynamicLogLog6.CALC_FGRA=Parse.parseBoolean(b);}
		else if(a.equals("saturate") || a.equals("sat")){UltraDynamicLogLog6.SATURATE_ON_OVERFLOW=Parse.parseBoolean(b);}
		else if(a.equals("printdlctiers")){DDLCalibrationDriver.PRINT_DLC_TIERS=Parse.parseBoolean(b);}
		else if(a.equals("frozenhistory") || a.equals("frozen")){UltraLogLog8.FROZEN_HISTORY=Parse.parseBoolean(b);}
		else if(a.equals("statepower") || a.equals("sp")){UltraLogLog8.STATE_POWER=Double.parseDouble(b);}
		else if(a.equals("statecfoffset") || a.equals("sco")){UltraLogLog8.STATE_CF_OFFSET=Double.parseDouble(b);}

		// PLL family flags
		else if(a.equals("plloffset") || a.equals("pco")){
			ProtoLogLog16b.CF_OFFSET=Double.parseDouble(b);
			ProtoLogLog16.CF_OFFSET=Double.parseDouble(b);
			ProtoLogLog16c.CF_OFFSET=Double.parseDouble(b);
		}
		else if(a.equals("hbits")){int v=Integer.parseInt(b); ProtoLogLog16.HISTORY_BITS=v; ProtoLogLog16b.HISTORY_BITS=v; ProtoLogLog16c.HISTORY_BITS=v;}
		else if(a.equals("lbits")){int v=Integer.parseInt(b); ProtoLogLog16.LUCK_BITS=v; ProtoLogLog16b.LUCK_BITS=v; ProtoLogLog16c.LUCK_BITS=v;}
		else if(a.equals("mbits")){int v=Integer.parseInt(b); ProtoLogLog16.MANTISSA_BITS=v; ProtoLogLog16b.MANTISSA_BITS=v; ProtoLogLog16c.MANTISSA_BITS=v;}
		else if(a.equals("abits")){int v=Integer.parseInt(b); ProtoLogLog16.ANDTISSA_BITS=v; ProtoLogLog16b.ANDTISSA_BITS=v; ProtoLogLog16c.ANDTISSA_BITS=v;}
		else if(a.equals("nbits")){int v=Integer.parseInt(b); ProtoLogLog16.NLZ2_BITS=v; ProtoLogLog16b.NLZ2_BITS=v; ProtoLogLog16c.NLZ2_BITS=v;}
		else if(a.equals("emode") || a.equals("estimatemode")){ProtoLogLog16b.ESTIMATE_MODE=Integer.parseInt(b);}
		else if(a.equals("nocf")){if(Parse.parseBoolean(b)){ProtoLogLog16b.CORRECTION_TABLE=null; ProtoLogLog16.CORRECTION_TABLE=null;}}

		// DynamicDemiLog8 flags
		else if(a.equals("empiricalmantissa") || a.equals("em")){DynamicDemiLog8.USE_EMPIRICAL_MANTISSA=Parse.parseBoolean(b);}
		else if(a.equals("mantissacfoffset") || a.equals("mco")){DynamicDemiLog8.MANTISSA_CF_OFFSET=Double.parseDouble(b);}

		// FLL2 flags
		else if(a.equals("fll2mult")){FutureLogLog2.TERMINAL_CORRECTION=Double.parseDouble(b);}
		else if(a.equals("clampoverflow")){FutureLogLog2.CLAMP_OVERFLOW=Parse.parseBoolean(b);}
		else if(a.equals("hcweight") || a.equals("ldlcweight")){CardinalityTracker.LDLC_HC_WEIGHT=Double.parseDouble(b);}

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
			DynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
			BankedDynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
			DualHashDynamicLogLog3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
			DynamicLogLog3v3.CORRECT_OVERFLOW=Parse.parseBoolean(b);
			DynamicLogLog3v4.CORRECT_OVERFLOW=Parse.parseBoolean(b);
			if(Parse.parseBoolean(b)){DynamicLogLog3.IGNORE_OVERFLOW=false; DynamicLogLog2.IGNORE_OVERFLOW=false; BankedDynamicLogLog3.IGNORE_OVERFLOW=false; DynamicLogLog3v2.IGNORE_OVERFLOW=false; DynamicLogLog3v3.IGNORE_OVERFLOW=false; DynamicLogLog3v4.IGNORE_OVERFLOW=false;}
		}
		else if(a.equals("ignoreoverflow") || a.equals("io")){
			DynamicLogLog3.IGNORE_OVERFLOW=Parse.parseBoolean(b);
			DynamicLogLog2.IGNORE_OVERFLOW=Parse.parseBoolean(b);
			BankedDynamicLogLog3.IGNORE_OVERFLOW=Parse.parseBoolean(b);
			DualHashDynamicLogLog3.IGNORE_OVERFLOW=Parse.parseBoolean(b);
			DynamicLogLog3v2.IGNORE_OVERFLOW=Parse.parseBoolean(b);
			DynamicLogLog3v3.IGNORE_OVERFLOW=Parse.parseBoolean(b);
			DynamicLogLog3v4.IGNORE_OVERFLOW=Parse.parseBoolean(b);
			if(Parse.parseBoolean(b)){DynamicLogLog3.CORRECT_OVERFLOW=false; BankedDynamicLogLog3.CORRECT_OVERFLOW=false; DynamicLogLog3v3.CORRECT_OVERFLOW=false; DynamicLogLog3v4.CORRECT_OVERFLOW=false;}
		}
		else if(a.equals("iobias") || a.equals("iob")){
			DynamicLogLog3.USE_IO_BIAS=Parse.parseBoolean(b);
			DynamicLogLog3v2.USE_IO_BIAS=Parse.parseBoolean(b);
			DynamicLogLog3v3.USE_IO_BIAS=Parse.parseBoolean(b);
			DynamicLogLog3v4.USE_IO_BIAS=Parse.parseBoolean(b);
			BankedDynamicLogLog3.USE_IO_BIAS=Parse.parseBoolean(b);
		}
		else if(a.equals("iobfile") || a.equals("iobf")){
			DynamicLogLog3.loadIOBias(b);
			DynamicLogLog3v2.loadIOBias(b);
			DynamicLogLog3v3.loadIOBias(b);
			DynamicLogLog3v4.loadIOBias(b);
			BankedDynamicLogLog3.loadIOBias(b);
		}
		else if(a.equals("clamptoadded") || a.equals("clamp")){CardinalityTracker.clampToAdded=Parse.parseBoolean(b);}
		else if(a.equals("dual")){
			DualHashDynamicLogLog3.DUAL=Parse.parseBoolean(b);
		}
		else{return false;}
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
	 * @param hcWeightExplicit true if hcweight was explicitly set (skip auto-defaults)
	 */
	public static void initializeAll(String loglogtype, int buckets, int k,
			String cffile, String pllmode, boolean hcWeightExplicit){

		// PLL mode application
		if(pllmode!=null){
			if(pllmode.equals("mantissa")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_MANTISSA); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_MANTISSA); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_MANTISSA);}
			else if(pllmode.equals("andtissa")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_ANDTISSA); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_ANDTISSA); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_ANDTISSA);}
			else if(pllmode.equals("nlz2")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_NLZ2); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_NLZ2); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_NLZ2);}
			else if(pllmode.equals("history")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_HISTORY); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_HISTORY); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_HISTORY);}
			else if(pllmode.equals("luck")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_LUCK); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_LUCK); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_LUCK);}
			else if(pllmode.equals("none")){ProtoLogLog16.setMode(ProtoLogLog16.MODE_NONE); ProtoLogLog16b.setMode(ProtoLogLog16b.MODE_NONE); ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_NONE);}
			else if(pllmode.equals("histmant") || pllmode.equals("historymantissa")){ProtoLogLog16c.setMode(ProtoLogLog16c.MODE_HISTORY|ProtoLogLog16c.MODE_MANTISSA);}
			else{throw new RuntimeException("Unknown pllmode: "+pllmode);}
		}

		// HC weight auto-defaults based on history bits
		if(!hcWeightExplicit){
			final int hb=ProtoLogLog16c.HISTORY_BITS;
			if(hb==1){CardinalityTracker.LDLC_HC_WEIGHT=0.456;}
			else if(hb==2){CardinalityTracker.LDLC_HC_WEIGHT=0.50;}
			else if(hb==3){CardinalityTracker.LDLC_HC_WEIGHT=0.475;}
		}

		// Class init and column whitelist
		DDLCalibrationDriver.makeInstance(loglogtype, buckets, k, 0L, 0);
		DDLCalibrationDriver.v3ColsForType(loglogtype);

		// Per-class Mean CF formula coefficients
		boolean classHasMeanCf=true, classHasHcCf=false, classHasHmeanmCf=false;
		if(loglogtype.equals("cdll4")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=1.5;}
		else if(loglogtype.equals("dhdll3")){CorrectionFactor.meanCfCoeffs=CorrectionFactor.MCF_DLL4; AbstractCardStats.TIER_SCALE=(DualHashDynamicLogLog3.DUAL ? 2 : 1.5);}
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
		else if(loglogtype.equals("udll6")){
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
		if("dll4".equals(loglogtype) || "dll4m".equals(loglogtype)){DynamicLogLog4.loadWordTable();}

		// Load CF table: explicit cffile overrides, otherwise auto-select per type+mode
		if(cffile==null){cffile=defaultCFFile(loglogtype);}
		if(cffile!=null){
			CorrectionFactor.initialize(cffile, buckets);
			if("pll16b".equals(loglogtype)){ProtoLogLog16b.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("pll16c".equals(loglogtype)){ProtoLogLog16c.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("udll6".equals(loglogtype)){UltraDynamicLogLog6.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
			else if("ertl".equals(loglogtype)){ErtlULL.setCFMatrix(CorrectionFactor.CF_MATRIX, buckets);}
		}

		// Publish immutable CF snapshot for worker threads
		CorrectionFactor.publishSnapshot();
	}

	/**
	 * Select the default CF table file for a given estimator type and overflow mode.
	 * Returns null if no table exists for the type.
	 */
	private static String defaultCFFile(String loglogtype){
		if(loglogtype.equals("cdll4") || loglogtype.equals("dhdll3") || loglogtype.equals("dhdll4")){
			return null; // No CF table yet for compressed/dual-hash variants
		}else if(loglogtype.equals("dll3")){
			if(DynamicLogLog3.IGNORE_OVERFLOW){return "?cardinalityCorrectionDLL3_iot.tsv.gz";}
			else if(!DynamicLogLog3.CORRECT_OVERFLOW){return "?cardinalityCorrectionDLL3_iof.tsv.gz";}
			else{return DynamicLogLog3.CF_FILE;}
		}else if(loglogtype.equals("dll4") || loglogtype.equals("dll4m")){
			return DynamicLogLog4.CF_FILE;
		}else if(loglogtype.equals("dll2")){
			if(DynamicLogLog2.IGNORE_OVERFLOW){return "?cardinalityCorrectionDLL2_iot.tsv.gz";}
			else{return "?cardinalityCorrectionDLL2_iof.tsv.gz";}
		}else if(loglogtype.equals("dll3v4")){
			if(DynamicLogLog3v4.IGNORE_OVERFLOW){return "?cardinalityCorrectionDLL3_iot.tsv.gz";}
			else if(!DynamicLogLog3v4.CORRECT_OVERFLOW){return "?cardinalityCorrectionDLL3_iof.tsv.gz";}
			else{return DynamicLogLog3v4.CF_FILE;}
		}else if(loglogtype.equals("dll3v2") || loglogtype.equals("dll3v3")){
			return DynamicLogLog3.CF_FILE;
		}else if(loglogtype.equals("bdll3")){
			if(BankedDynamicLogLog3.IGNORE_OVERFLOW){return "?cardinalityCorrectionBDLL3_iot.tsv.gz";}
			else if(!BankedDynamicLogLog3.CORRECT_OVERFLOW){return "?cardinalityCorrectionBDLL3_cof.tsv.gz";}
			else{return "?cardinalityCorrectionBDLL3_cot.tsv.gz";}
		}else if(loglogtype.equals("ll6")){
			return LogLog6.CF_FILE;
		}else if(loglogtype.equals("udll6")){
			return UltraDynamicLogLog6.CF_FILE;
		}else if(loglogtype.equals("ddl") || loglogtype.equals("ddl10")){
			return DynamicDemiLog.CF_FILE;
		}else if(loglogtype.equals("ddl8")){
			return DynamicDemiLog8.CF_FILE;
		}else if(loglogtype.equals("ttll")){
			return "?"+TwinTailLogLog.CF_FILE;
		}
		return null;
	}

}
