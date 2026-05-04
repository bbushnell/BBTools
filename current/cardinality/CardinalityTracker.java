package cardinality;

import java.util.Arrays;
import shared.Random;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import jgi.Dedupe;
import parse.Parser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import structures.LongList;
import structures.SuperLongList;
import ukmer.Kmer;

/**
 * Abstract superclass for cardinality-tracking structures like LogLog.
 * Provides probabilistic estimation of unique k-mer counts in sequences using various LogLog-based algorithms.
 * Supports different implementations optimized for accuracy, speed, or memory usage.
 * @author Brian Bushnell
 * @date Feb 20, 2020
 */
public abstract class CardinalityTracker implements Drivable {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private static final String pickType(final String type){
		if(trackCounts){
			if("BBLog".equalsIgnoreCase(type)){return type;}
			else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type) ||
				"DDL".equalsIgnoreCase(type) || "DynamicDemiLog".equalsIgnoreCase(type)){
				return type;
			}else{return "LogLog16";}
		}
		return type;
	}

	/**
	 * Factory method that creates a tracker using default settings.
	 * Subclass is determined by static Parser.loglogType field.
	 * BBLog is preferred when trackCounts is enabled for optimal accuracy and speed.
	 * @return New CardinalityTracker instance of the configured type
	 */
	public static CardinalityTracker makeTracker(final String type0){
		assert(type0==null || type0.equalsIgnoreCase(Parser.loglogType)) :
			"Passed type '"+type0+"' differs from configured Parser.loglogType '"+Parser.loglogType+"'";
		final String type=pickType(type0);
		if("BBLog".equalsIgnoreCase(type)){
			return new BBLog();//Fastest, most accurate
		}else if("LogLog".equalsIgnoreCase(type)){
			return new LogLog();//Least accurate
		}else if("LogLog2".equalsIgnoreCase(type)){
			return new LogLog2();//Slowest, uses mantissa
		}else if("LogLog16".equalsIgnoreCase(type) || "LL16".equalsIgnoreCase(type)){
			return new LogLog16();//Uses 10-bit mantissa
//		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type)){
//			return new DynamicLogLog();//Uses 10-bit mantissa
		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type) || 
			"DDL".equalsIgnoreCase(type) || "DynamicDemiLog".equalsIgnoreCase(type)){
			return new DynamicDemiLog();//Uses 10-bit mantissa
		}else if("LogLog8".equalsIgnoreCase(type)){
			return new LogLog8();//Lowest memory
		}else if("InvertedLogLog".equalsIgnoreCase(type) || "ILL".equalsIgnoreCase(type)){
			return new InvertedLogLog();//Experimental: stores max low-bits per NLZ bucket
		}else if("DDL2".equalsIgnoreCase(type) || "DynamicDemiLog2".equalsIgnoreCase(type)){
			return new DynamicDemiLog2();
		}else if("DDL8".equalsIgnoreCase(type) || "DynamicDemiLog8".equalsIgnoreCase(type)){
			return new DynamicDemiLog8();
		}else if("DDL8v2".equalsIgnoreCase(type) || "DynamicDemiLog8v2".equalsIgnoreCase(type)){
			return new DynamicDemiLog8v2();
		}else if("CDLL4".equalsIgnoreCase(type) || "CompressedDynamicLogLog4".equalsIgnoreCase(type)){
			return new CompressedDynamicLogLog4();
		}else if("VCDLL4".equalsIgnoreCase(type) || "VariableCompressedDynamicLogLog4".equalsIgnoreCase(type)){
			return new VariableCompressedDynamicLogLog4();
		}else if("CDLL5".equalsIgnoreCase(type) || "CompressedDynamicLogLog5".equalsIgnoreCase(type)){
			return new CompressedDynamicLogLog5();
		}else if("BCDLL5".equalsIgnoreCase(type) || "BankedCompressedDynamicLogLog5".equalsIgnoreCase(type)){
			return new BankedCompressedDynamicLogLog5();
		}else if("ACDLL5".equalsIgnoreCase(type) || "ArithmeticCompressedDynamicLogLog5".equalsIgnoreCase(type)){
			return new ArithmeticCompressedDynamicLogLog5();
		}else if("ACDLL4".equalsIgnoreCase(type) || "ArithmeticCompressedDynamicLogLog4".equalsIgnoreCase(type)){
			return new ArithmeticCompressedDynamicLogLog4();
		}else if("HCDLL4".equalsIgnoreCase(type) || "HalfCappedDynamicLogLog4".equalsIgnoreCase(type)){
			return new HalfCappedDynamicLogLog4();
		}else if("EDLL8".equalsIgnoreCase(type) || "ExpandedDynamicLogLog8".equalsIgnoreCase(type)){
			return new ExpandedDynamicLogLog8();
		}else if("EDLL9".equalsIgnoreCase(type) || "ExpandedDynamicLogLog9".equalsIgnoreCase(type)){
			return new ExpandedDynamicLogLog9();
		}else if("EXA".equalsIgnoreCase(type) || "ExaLogLog".equalsIgnoreCase(type)){
			return makeByReflection("ExaLogLogWrapper");
		}else if("AUDLL32".equalsIgnoreCase(type) || "ArithmeticUltraDynamicLogLog32".equalsIgnoreCase(type)){
			return new ArithmeticUltraDynamicLogLog32();
		}else if("AUDLL33".equalsIgnoreCase(type) || "ArithmeticUltraDynamicLogLog33".equalsIgnoreCase(type)){
			return new ArithmeticUltraDynamicLogLog33();
		}else if("AVDLL32".equalsIgnoreCase(type) || "ArithmeticVariableDynamicLogLog32".equalsIgnoreCase(type)){
			return new ArithmeticVariableDynamicLogLog32();
		}else if("AVDLL36".equalsIgnoreCase(type) || "ArithmeticVariableDynamicLogLog36".equalsIgnoreCase(type)){
			return new ArithmeticVariableDynamicLogLog36();
		}else if("AVDLL34".equalsIgnoreCase(type) || "ArithmeticVariableDynamicLogLog34".equalsIgnoreCase(type)){
			return new ArithmeticVariableDynamicLogLog34();
		}else if("DHDLL3".equalsIgnoreCase(type) || "DualHashDynamicLogLog3".equalsIgnoreCase(type) || "cdll3".equalsIgnoreCase(type) || "CompressedDynamicLogLog3".equalsIgnoreCase(type)){
			return new CompressedDynamicLogLog3();
		}else if("DHDLL4".equalsIgnoreCase(type) || "DualHashDynamicLogLog4".equalsIgnoreCase(type)){
			return new DualHashDynamicLogLog4();
		}else if("DLL4".equalsIgnoreCase(type) || "DynamicLogLog4".equalsIgnoreCase(type) ||
			"DDL4".equalsIgnoreCase(type) || "DynamicDemiLog4".equalsIgnoreCase(type)){
			return new DynamicLogLog4();
		}else if("DLL3".equalsIgnoreCase(type) || "DynamicLogLog3".equalsIgnoreCase(type) ||
			"DDL3".equalsIgnoreCase(type) || "DynamicDemiLog3".equalsIgnoreCase(type)){
			return new DynamicLogLog3();
		}else if("DLL3v4".equalsIgnoreCase(type) || "DynamicLogLog3v4".equalsIgnoreCase(type)){
			return new DynamicLogLog3v4();
		}else if("DLL2".equalsIgnoreCase(type) || "DynamicLogLog2".equalsIgnoreCase(type)){
			return new DynamicLogLog2();
		}else if("LL6".equalsIgnoreCase(type) || "LogLog6".equalsIgnoreCase(type)){
			return new LogLog6();
		}else if("ULL".equalsIgnoreCase(type) || "ErtlULL".equalsIgnoreCase(type) || "UltraLogLog".equalsIgnoreCase(type)){
			return new ErtlULL();
		}else if("UDLL6".equalsIgnoreCase(type) || "UltraDynamicLogLog6".equalsIgnoreCase(type)){
			return new UltraDynamicLogLog6();
		}else if("UDLL36".equalsIgnoreCase(type) || "UltraDynamicLogLog36".equalsIgnoreCase(type)){
			return new UltraDynamicLogLog36();
		}else if("BDLL3".equalsIgnoreCase(type) || "BankedDynamicLogLog3".equalsIgnoreCase(type)){
			return new BankedDynamicLogLog3();
		}else if("BCDLL3".equalsIgnoreCase(type) || "BankedCompressedDynamicLogLog3".equalsIgnoreCase(type)){
			return new BankedCompressedDynamicLogLog3();
		}else if("BDLL4".equalsIgnoreCase(type) || "BankedDynamicLogLog4".equalsIgnoreCase(type)){
			return new BankedDynamicLogLog4();
		}else if("BDLL5".equalsIgnoreCase(type) || "BankedDynamicLogLog5".equalsIgnoreCase(type)){
			return new BankedDynamicLogLog5();
		}else if("PLL16c".equalsIgnoreCase(type)){
			return new ProtoLogLog16c();
		}else if("ErtlB".equalsIgnoreCase(type) || "ErtlULLb".equalsIgnoreCase(type)){
			return new ErtlULLb();
		}else if("DLL4m".equalsIgnoreCase(type) || "DynamicLogLog4m".equalsIgnoreCase(type)){
			return new DynamicLogLog4m();
		}else if("HTB".equalsIgnoreCase(type) || "HyperTwoBits".equalsIgnoreCase(type)){
			return new HyperTwoBits();
		}else if("HTC".equalsIgnoreCase(type) || "HLLTailCut".equalsIgnoreCase(type)){
			return new HLLTailCut();
		}else if("HTC4".equalsIgnoreCase(type) || "HLLTailCut4".equalsIgnoreCase(type)){
			return new HLLTailCut4();
		}else if("FLL2".equalsIgnoreCase(type) || "FutureLogLog2".equalsIgnoreCase(type)){
			return new FutureLogLog2();
		}else if("BCLL".equalsIgnoreCase(type) || "BankedCeilingLogLog".equalsIgnoreCase(type)){
			return new BankedCeilingLogLog();
		}else if("TTLL".equalsIgnoreCase(type) || "TwinTailLogLog".equalsIgnoreCase(type)){
			return new TwinTailLogLog();
		}else if("TTLL3".equalsIgnoreCase(type) || "TwinTailLogLog3".equalsIgnoreCase(type)){
			return new TwinTailLogLog3();
		}else if("TTLL4".equalsIgnoreCase(type) || "TwinTailLogLog4".equalsIgnoreCase(type)){
			return new TwinTailLogLog4();
		}else if("CTTLL".equalsIgnoreCase(type) || "CompressedTwinTailLogLog".equalsIgnoreCase(type)){
			return new CompressedTwinTailLogLog();
		}
		assert(false) : "TODO: "+type;
		throw new RuntimeException(type);
	}

	/**
	 * Factory method that creates a tracker using default settings.
	 * Subclass is determined by static Parser.loglogType field.
	 * BBLog is preferred when trackCounts is enabled for optimal accuracy and speed.
	 * @return New CardinalityTracker instance of the configured type
	 */
	public static CardinalityTracker makeTracker(){return makeTracker(Parser.loglogType);}

	/**
	 * Factory method that creates a tracker using parsed settings.
	 * Subclass is determined by static type field.
	 * Parameters are extracted from the Parser object.
	 * @param p Parser containing configuration parameters
	 * @return New CardinalityTracker instance configured from parser
	 */
	public static CardinalityTracker makeTracker(final Parser p){
		final String type=pickType(Parser.loglogType);
		if("BBLog".equalsIgnoreCase(type)){
			return new BBLog(p);
		}else if("LogLog".equalsIgnoreCase(type)){
			return new LogLog(p);
		}else if("LogLog2".equalsIgnoreCase(type)){
			return new LogLog2(p);
		}else if("LogLog16".equalsIgnoreCase(type) || "LL16".equalsIgnoreCase(type)){
			return new LogLog16(p);
//		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type)){
//			return new DynamicLogLog(p);
		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type) || 
			"DDL".equalsIgnoreCase(type) || "DynamicDemiLog".equalsIgnoreCase(type)){
			return new DynamicDemiLog(p);
		}else if("LogLog8".equalsIgnoreCase(type)){
			return new LogLog8(p);
		}else if("InvertedLogLog".equalsIgnoreCase(type) || "ILL".equalsIgnoreCase(type)){
			return new InvertedLogLog(p);//Experimental: stores max low-bits per NLZ bucket
		}else if("DDL2".equalsIgnoreCase(type) || "DynamicDemiLog2".equalsIgnoreCase(type)){
			return new DynamicDemiLog2(p);
		}else if("DDL8".equalsIgnoreCase(type) || "DynamicDemiLog8".equalsIgnoreCase(type)){
			return new DynamicDemiLog8(p);
		}else if("DDL8v2".equalsIgnoreCase(type) || "DynamicDemiLog8v2".equalsIgnoreCase(type)){
			return new DynamicDemiLog8v2(p);
		}else if("CDLL4".equalsIgnoreCase(type) || "CompressedDynamicLogLog4".equalsIgnoreCase(type)){
			return new CompressedDynamicLogLog4(p);
		}else if("VCDLL4".equalsIgnoreCase(type) || "VariableCompressedDynamicLogLog4".equalsIgnoreCase(type)){
			return new VariableCompressedDynamicLogLog4(p);
		}else if("CDLL5".equalsIgnoreCase(type) || "CompressedDynamicLogLog5".equalsIgnoreCase(type)){
			return new CompressedDynamicLogLog5(p);
		}else if("BCDLL5".equalsIgnoreCase(type) || "BankedCompressedDynamicLogLog5".equalsIgnoreCase(type)){
			return new BankedCompressedDynamicLogLog5(p);
		}else if("ACDLL5".equalsIgnoreCase(type) || "ArithmeticCompressedDynamicLogLog5".equalsIgnoreCase(type)){
			return new ArithmeticCompressedDynamicLogLog5(p);
		}else if("ACDLL4".equalsIgnoreCase(type) || "ArithmeticCompressedDynamicLogLog4".equalsIgnoreCase(type)){
			return new ArithmeticCompressedDynamicLogLog4(p);
		}else if("HCDLL4".equalsIgnoreCase(type) || "HalfCappedDynamicLogLog4".equalsIgnoreCase(type)){
			return new HalfCappedDynamicLogLog4(p);
		}else if("EDLL8".equalsIgnoreCase(type) || "ExpandedDynamicLogLog8".equalsIgnoreCase(type)){
			return new ExpandedDynamicLogLog8(p);
		}else if("EDLL9".equalsIgnoreCase(type) || "ExpandedDynamicLogLog9".equalsIgnoreCase(type)){
			return new ExpandedDynamicLogLog9(p);
		}else if("EXA".equalsIgnoreCase(type) || "ExaLogLog".equalsIgnoreCase(type)){
			return makeByReflection("ExaLogLogWrapper", p);
		}else if("AUDLL32".equalsIgnoreCase(type) || "ArithmeticUltraDynamicLogLog32".equalsIgnoreCase(type)){
			return new ArithmeticUltraDynamicLogLog32(p);
		}else if("AUDLL33".equalsIgnoreCase(type) || "ArithmeticUltraDynamicLogLog33".equalsIgnoreCase(type)){
			return new ArithmeticUltraDynamicLogLog33(p);
		}else if("AVDLL32".equalsIgnoreCase(type) || "ArithmeticVariableDynamicLogLog32".equalsIgnoreCase(type)){
			return new ArithmeticVariableDynamicLogLog32(p);
		}else if("AVDLL36".equalsIgnoreCase(type) || "ArithmeticVariableDynamicLogLog36".equalsIgnoreCase(type)){
			return new ArithmeticVariableDynamicLogLog36(p);
		}else if("AVDLL34".equalsIgnoreCase(type) || "ArithmeticVariableDynamicLogLog34".equalsIgnoreCase(type)){
			return new ArithmeticVariableDynamicLogLog34(p);
		}else if("DHDLL3".equalsIgnoreCase(type) || "DualHashDynamicLogLog3".equalsIgnoreCase(type) || "cdll3".equalsIgnoreCase(type) || "CompressedDynamicLogLog3".equalsIgnoreCase(type)){
			return new CompressedDynamicLogLog3(p);
		}else if("DHDLL4".equalsIgnoreCase(type) || "DualHashDynamicLogLog4".equalsIgnoreCase(type)){
			return new DualHashDynamicLogLog4(p);
		}else if("DLL4".equalsIgnoreCase(type) || "DynamicLogLog4".equalsIgnoreCase(type) ||
			"DDL4".equalsIgnoreCase(type) || "DynamicDemiLog4".equalsIgnoreCase(type)){
			return new DynamicLogLog4(p);
		}else if("DLL3".equalsIgnoreCase(type) || "DynamicLogLog3".equalsIgnoreCase(type) ||
			"DDL3".equalsIgnoreCase(type) || "DynamicDemiLog3".equalsIgnoreCase(type)){
			return new DynamicLogLog3(p);
		}else if("DLL3v4".equalsIgnoreCase(type) || "DynamicLogLog3v4".equalsIgnoreCase(type)){
			return new DynamicLogLog3v4(p);
		}else if("DLL2".equalsIgnoreCase(type) || "DynamicLogLog2".equalsIgnoreCase(type)){
			return new DynamicLogLog2(p);
		}else if("LL6".equalsIgnoreCase(type) || "LogLog6".equalsIgnoreCase(type)){
			return new LogLog6(p);
		}else if("ULL".equalsIgnoreCase(type) || "ErtlULL".equalsIgnoreCase(type) || "UltraLogLog".equalsIgnoreCase(type)){
			return new ErtlULL(p);
		}else if("UDLL6".equalsIgnoreCase(type) || "UltraDynamicLogLog6".equalsIgnoreCase(type)){
			return new UltraDynamicLogLog6(p);
		}else if("UDLL36".equalsIgnoreCase(type) || "UltraDynamicLogLog36".equalsIgnoreCase(type)){
			return new UltraDynamicLogLog36(p);
		}else if("BDLL3".equalsIgnoreCase(type) || "BankedDynamicLogLog3".equalsIgnoreCase(type)){
			return new BankedDynamicLogLog3(p);
		}else if("BCDLL3".equalsIgnoreCase(type) || "BankedCompressedDynamicLogLog3".equalsIgnoreCase(type)){
			return new BankedCompressedDynamicLogLog3(p);
		}else if("BDLL4".equalsIgnoreCase(type) || "BankedDynamicLogLog4".equalsIgnoreCase(type)){
			return new BankedDynamicLogLog4(p);
		}else if("BDLL5".equalsIgnoreCase(type) || "BankedDynamicLogLog5".equalsIgnoreCase(type)){
			return new BankedDynamicLogLog5(p);
		}else if("PLL16c".equalsIgnoreCase(type)){
			return new ProtoLogLog16c(p);
		}else if("ErtlB".equalsIgnoreCase(type) || "ErtlULLb".equalsIgnoreCase(type)){
			return new ErtlULLb(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
		}else if("DLL4m".equalsIgnoreCase(type) || "DynamicLogLog4m".equalsIgnoreCase(type)){
			return new DynamicLogLog4m(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
		}else if("HTB".equalsIgnoreCase(type) || "HyperTwoBits".equalsIgnoreCase(type)){
			return new HyperTwoBits(p);
		}else if("HTC".equalsIgnoreCase(type) || "HLLTailCut".equalsIgnoreCase(type)){
			return new HLLTailCut(p);
		}else if("HTC4".equalsIgnoreCase(type) || "HLLTailCut4".equalsIgnoreCase(type)){
			return new HLLTailCut4(p);
		}else if("FLL2".equalsIgnoreCase(type) || "FutureLogLog2".equalsIgnoreCase(type)){
			return new FutureLogLog2(p);
		}else if("BCLL".equalsIgnoreCase(type) || "BankedCeilingLogLog".equalsIgnoreCase(type)){
			return new BankedCeilingLogLog(p);
		}else if("TTLL".equalsIgnoreCase(type) || "TwinTailLogLog".equalsIgnoreCase(type)){
			return new TwinTailLogLog(p);
		}else if("TTLL3".equalsIgnoreCase(type) || "TwinTailLogLog3".equalsIgnoreCase(type)){
			return new TwinTailLogLog3(p);
		}else if("TTLL4".equalsIgnoreCase(type) || "TwinTailLogLog4".equalsIgnoreCase(type)){
			return new TwinTailLogLog4(p);
		}else if("CTTLL".equalsIgnoreCase(type) || "CompressedTwinTailLogLog".equalsIgnoreCase(type)){
			return new CompressedTwinTailLogLog(p);
		}
		assert(false) : "TODO: "+type;
		throw new RuntimeException(type);
	}

	/**
	 * Factory method that creates a tracker with specified settings.
	 * Subclass is determined by static type field.
	 * Allows direct specification of all key parameters.
	 * @param buckets_ Number of buckets (will be rounded to next power of 2)
	 * @param k_ K-mer length for hashing
	 * @param seed Random number generator seed (-1 for random seed)
	 * @param minProb_ Ignore k-mers with correctness probability below this threshold
	 * @return New CardinalityTracker instance with specified configuration
	 */
	public static CardinalityTracker makeTracker(final int buckets_, final int k_, final long seed, final float minProb_){
		return makeTracker(Parser.loglogType, buckets_, k_, seed, minProb_);
	}

	/** Thread-safe factory: takes type explicitly instead of reading Parser.loglogType. */
	public static CardinalityTracker makeTracker(final String type0, final int buckets_, final int k_, final long seed, final float minProb_){
		final String type=pickType(type0);
		if("BBLog".equalsIgnoreCase(type)){
			return new BBLog(buckets_, k_, seed, minProb_);
		}else if("LogLog".equalsIgnoreCase(type)){
			return new LogLog(buckets_, k_, seed, minProb_);
		}else if("LogLog2".equalsIgnoreCase(type)){
			return new LogLog2(buckets_, k_, seed, minProb_);
		}else if("LogLog16".equalsIgnoreCase(type) || "LL16".equalsIgnoreCase(type)){
			return new LogLog16(buckets_, k_, seed, minProb_);
//		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type)) {
//			return new DynamicLogLog(buckets_, k_, seed, minProb_);
		}else if("DLL".equalsIgnoreCase(type) || "DynamicLogLog".equalsIgnoreCase(type) || 
			"DDL".equalsIgnoreCase(type) || "DynamicDemiLog".equalsIgnoreCase(type)){
			return new DynamicDemiLog(buckets_, k_, seed, minProb_);
		}else if("LogLog8".equalsIgnoreCase(type)){
			return new LogLog8(buckets_, k_, seed, minProb_);
		}else if("InvertedLogLog".equalsIgnoreCase(type) || "ILL".equalsIgnoreCase(type)){
			return new InvertedLogLog(buckets_, k_, seed, minProb_);//Experimental: stores max low-bits per NLZ bucket
		}else if("DDL2".equalsIgnoreCase(type) || "DynamicDemiLog2".equalsIgnoreCase(type)){
			return new DynamicDemiLog2(buckets_, k_, seed, minProb_);
		}else if("DDL8".equalsIgnoreCase(type) || "DynamicDemiLog8".equalsIgnoreCase(type)){
			return new DynamicDemiLog8(buckets_, k_, seed, minProb_);
		}else if("DDL8v2".equalsIgnoreCase(type) || "DynamicDemiLog8v2".equalsIgnoreCase(type)){
			return new DynamicDemiLog8v2(buckets_, k_, seed, minProb_);
		}else if("CDLL4".equalsIgnoreCase(type) || "CompressedDynamicLogLog4".equalsIgnoreCase(type)){
			return new CompressedDynamicLogLog4(buckets_, k_, seed, minProb_);
		}else if("VCDLL4".equalsIgnoreCase(type) || "VariableCompressedDynamicLogLog4".equalsIgnoreCase(type)){
			return new VariableCompressedDynamicLogLog4(buckets_, k_, seed, minProb_);
		}else if("CDLL5".equalsIgnoreCase(type) || "CompressedDynamicLogLog5".equalsIgnoreCase(type)){
			return new CompressedDynamicLogLog5(buckets_, k_, seed, minProb_);
		}else if("BCDLL5".equalsIgnoreCase(type) || "BankedCompressedDynamicLogLog5".equalsIgnoreCase(type)){
			return new BankedCompressedDynamicLogLog5(buckets_, k_, seed, minProb_);
		}else if("ACDLL5".equalsIgnoreCase(type) || "ArithmeticCompressedDynamicLogLog5".equalsIgnoreCase(type)){
			return new ArithmeticCompressedDynamicLogLog5(buckets_, k_, seed, minProb_);
		}else if("ACDLL4".equalsIgnoreCase(type) || "ArithmeticCompressedDynamicLogLog4".equalsIgnoreCase(type)){
			return new ArithmeticCompressedDynamicLogLog4(buckets_, k_, seed, minProb_);
		}else if("HCDLL4".equalsIgnoreCase(type) || "HalfCappedDynamicLogLog4".equalsIgnoreCase(type)){
			return new HalfCappedDynamicLogLog4(buckets_, k_, seed, minProb_);
		}else if("EDLL8".equalsIgnoreCase(type) || "ExpandedDynamicLogLog8".equalsIgnoreCase(type)){
			return new ExpandedDynamicLogLog8(buckets_, k_, seed, minProb_);
		}else if("EDLL9".equalsIgnoreCase(type) || "ExpandedDynamicLogLog9".equalsIgnoreCase(type)){
			return new ExpandedDynamicLogLog9(buckets_, k_, seed, minProb_);
		}else if("EXA".equalsIgnoreCase(type) || "ExaLogLog".equalsIgnoreCase(type)){
			return makeByReflection("ExaLogLogWrapper", buckets_, k_, seed, minProb_);
		}else if("AUDLL32".equalsIgnoreCase(type) || "ArithmeticUltraDynamicLogLog32".equalsIgnoreCase(type)){
			return new ArithmeticUltraDynamicLogLog32(buckets_, k_, seed, minProb_);
		}else if("AUDLL33".equalsIgnoreCase(type) || "ArithmeticUltraDynamicLogLog33".equalsIgnoreCase(type)){
			return new ArithmeticUltraDynamicLogLog33(buckets_, k_, seed, minProb_);
		}else if("AVDLL32".equalsIgnoreCase(type) || "ArithmeticVariableDynamicLogLog32".equalsIgnoreCase(type)){
			return new ArithmeticVariableDynamicLogLog32(buckets_, k_, seed, minProb_);
		}else if("AVDLL36".equalsIgnoreCase(type) || "ArithmeticVariableDynamicLogLog36".equalsIgnoreCase(type)){
			return new ArithmeticVariableDynamicLogLog36(buckets_, k_, seed, minProb_);
		}else if("AVDLL34".equalsIgnoreCase(type) || "ArithmeticVariableDynamicLogLog34".equalsIgnoreCase(type)){
			return new ArithmeticVariableDynamicLogLog34(buckets_, k_, seed, minProb_);
		}else if("DHDLL3".equalsIgnoreCase(type) || "DualHashDynamicLogLog3".equalsIgnoreCase(type) || "cdll3".equalsIgnoreCase(type) || "CompressedDynamicLogLog3".equalsIgnoreCase(type)){
			return new CompressedDynamicLogLog3(buckets_, k_, seed, minProb_);
		}else if("DHDLL4".equalsIgnoreCase(type) || "DualHashDynamicLogLog4".equalsIgnoreCase(type)){
			return new DualHashDynamicLogLog4(buckets_, k_, seed, minProb_);
		}else if("DLL4".equalsIgnoreCase(type) || "DynamicLogLog4".equalsIgnoreCase(type) ||
			"DDL4".equalsIgnoreCase(type) || "DynamicDemiLog4".equalsIgnoreCase(type)){
			return new DynamicLogLog4(buckets_, k_, seed, minProb_);
		}else if("DLL3".equalsIgnoreCase(type) || "DynamicLogLog3".equalsIgnoreCase(type) ||
			"DDL3".equalsIgnoreCase(type) || "DynamicDemiLog3".equalsIgnoreCase(type)){
			return new DynamicLogLog3(buckets_, k_, seed, minProb_);
		}else if("DLL3v4".equalsIgnoreCase(type) || "DynamicLogLog3v4".equalsIgnoreCase(type)){
			return new DynamicLogLog3v4(buckets_, k_, seed, minProb_);
		}else if("DLL2".equalsIgnoreCase(type) || "DynamicLogLog2".equalsIgnoreCase(type)){
			return new DynamicLogLog2(buckets_, k_, seed, minProb_);
		}else if("LL6".equalsIgnoreCase(type) || "LogLog6".equalsIgnoreCase(type)){
			return new LogLog6(buckets_, k_, seed, minProb_);
		}else if("ULL".equalsIgnoreCase(type) || "ErtlULL".equalsIgnoreCase(type) || "UltraLogLog".equalsIgnoreCase(type)){
			return new ErtlULL(buckets_, k_, seed, minProb_);
		}else if("UDLL6".equalsIgnoreCase(type) || "UltraDynamicLogLog6".equalsIgnoreCase(type)){
			return new UltraDynamicLogLog6(buckets_, k_, seed, minProb_);
		}else if("UDLL6m".equalsIgnoreCase(type) || "UltraDynamicLogLog6m".equalsIgnoreCase(type)){
			return new UltraDynamicLogLog6m(buckets_, k_, seed, minProb_);
		}else if("UDLL36".equalsIgnoreCase(type) || "UltraDynamicLogLog36".equalsIgnoreCase(type)){
			return new UltraDynamicLogLog36(buckets_, k_, seed, minProb_);
		}else if("BDLL3".equalsIgnoreCase(type) || "BankedDynamicLogLog3".equalsIgnoreCase(type)){
			return new BankedDynamicLogLog3(buckets_, k_, seed, minProb_);
		}else if("BCDLL3".equalsIgnoreCase(type) || "BankedCompressedDynamicLogLog3".equalsIgnoreCase(type)){
			return new BankedCompressedDynamicLogLog3(buckets_, k_, seed, minProb_);
		}else if("BDLL4".equalsIgnoreCase(type) || "BankedDynamicLogLog4".equalsIgnoreCase(type)){
			return new BankedDynamicLogLog4(buckets_, k_, seed, minProb_);
		}else if("BDLL5".equalsIgnoreCase(type) || "BankedDynamicLogLog5".equalsIgnoreCase(type)){
			return new BankedDynamicLogLog5(buckets_, k_, seed, minProb_);
		}else if("PLL16c".equalsIgnoreCase(type)){
			return new ProtoLogLog16c(buckets_, k_, seed, minProb_);
		}else if("ErtlB".equalsIgnoreCase(type) || "ErtlULLb".equalsIgnoreCase(type)){
			return new ErtlULLb(buckets_, k_, seed, minProb_);
		}else if("DLL4m".equalsIgnoreCase(type) || "DynamicLogLog4m".equalsIgnoreCase(type)){
			return new DynamicLogLog4m(buckets_, k_, seed, minProb_);
		}else if("HTB".equalsIgnoreCase(type) || "HyperTwoBits".equalsIgnoreCase(type)){
			return new HyperTwoBits(buckets_, k_, seed, minProb_);
		}else if("HTC".equalsIgnoreCase(type) || "HLLTailCut".equalsIgnoreCase(type)){
			return new HLLTailCut(buckets_, k_, seed, minProb_);
		}else if("HTC4".equalsIgnoreCase(type) || "HLLTailCut4".equalsIgnoreCase(type)){
			return new HLLTailCut4(buckets_, k_, seed, minProb_);
		}else if("FLL2".equalsIgnoreCase(type) || "FutureLogLog2".equalsIgnoreCase(type)){
			return new FutureLogLog2(buckets_, k_, seed, minProb_);
		}else if("BCLL".equalsIgnoreCase(type) || "BankedCeilingLogLog".equalsIgnoreCase(type)){
			return new BankedCeilingLogLog(buckets_, k_, seed, minProb_);
		}else if("TTLL".equalsIgnoreCase(type) || "TwinTailLogLog".equalsIgnoreCase(type)){
			return new TwinTailLogLog(buckets_, k_, seed, minProb_);
		}else if("TTLL3".equalsIgnoreCase(type) || "TwinTailLogLog3".equalsIgnoreCase(type)){
			return new TwinTailLogLog3(buckets_, k_, seed, minProb_);
		}else if("TTLL4".equalsIgnoreCase(type) || "TwinTailLogLog4".equalsIgnoreCase(type)){
			return new TwinTailLogLog4(buckets_, k_, seed, minProb_);
		}else if("CTTLL".equalsIgnoreCase(type) || "CompressedTwinTailLogLog".equalsIgnoreCase(type)){
			return new CompressedTwinTailLogLog(buckets_, k_, seed, minProb_);
		}else if("ETTLL5".equalsIgnoreCase(type) || "ExpandedTwinTailLogLog5".equalsIgnoreCase(type)){
			return new ExpandedTwinTailLogLog5(buckets_, k_, seed, minProb_);
		}
		assert(false) : "TODO: "+type;
		throw new RuntimeException(type);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a tracker with parameters extracted from a Parser.
	 * @param p Parser containing loglog configuration values */
	public CardinalityTracker(final Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	/**
	 * Creates a tracker with specified parameters.
	 * Buckets will be rounded up to the next power of 2 for efficient bit masking.
	 * Initializes hash function with random XOR value from the specified seed.
	 * @param buckets_ Number of buckets (counters) - will be made power of 2
	 * @param k_ K-mer length for sequence hashing
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore k-mers with under this probability of being correct
	 */
	public CardinalityTracker(final int buckets_, final int k_, final long seed, final float minProb_){
		buckets=powerOf2AtLeast(buckets_);
		assert(buckets>0 && Integer.bitCount(buckets)==1) : "Buckets must be a power of 2: "+buckets;
		bucketMask=buckets-1;
		bucketBits=Long.bitCount(bucketMask);
		k=Kmer.getKbig(k_);
		minProb=minProb_;
		hashXor=(seed==-1 ? defaultSeed : seed);
//		assert(false) : seed+", "+randy+", "+hashXor;
		//Xor is not used anymore, but could be, as long as each copy gets the same one.
	}

	/** Create an independent copy of this tracker. */
	public abstract CardinalityTracker copy();

	/**
	 * Returns the lowest power of 2 that is greater than or equal to target.
	 * Required because buckets must be a power of 2 for efficient bit masking.
	 * @param target The minimum value needed
	 * @return Smallest power of 2 >= target, capped at 0x40000000
	 */
	public static final int powerOf2AtLeast(final int target){
		if(target<1){return 1;}
		int ret=1;
		final int limit=Tools.min(target, 0x40000000);
		while(ret<limit){ret<<=1;}
		return ret;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Program entry point (deprecated).
	 * Redirects to LogLogWrapper for actual processing.
	 * @param args Command-line arguments
	 */
	public static final void main(final String[] args){
		final LogLogWrapper llw=new LogLogWrapper(args);

		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;

		llw.process();

		Read.VALIDATE_IN_CONSTRUCTOR=vic;
	}

	/**
	 * Hashes and adds a number to this tracker.
	 * This is the ONLY correct way to add elements externally.
	 * Handles added-count tracking for clampToAdded support.
	 * @param number The value to hash and track
	 */
	public final void add(final long number){
		hashAndStore(number);
		added++;
	}

	/**
	 * Hashes and tracks all k-mers from a Read and its mate.
	 * Processes both forward and mate sequences if they meet minimum length requirements.
	 * @param r The Read to process (may be null)
	 */
	public final void hash(final Read r){
		if(r==null){return;}
		if(r.length()>=k){hash(r.bases, r.quality);}
		if(r.mateLength()>=k){hash(r.mate.bases, r.mate.quality);}
	}

	/**
	 * Hashes and tracks all k-mers from a SamLine.
	 * Processes the sequence if it meets minimum length requirements.
	 * @param r The SamLine to process (may be null)
	 */
	public final void hash(final SamLine r){
		if(r==null || r.seq==null){return;}
		if(r.length()>=k){hash(r.seq, r.qual);}
	}

	/**
	 * Hashes and tracks all k-mers from a sequence with quality scores.
	 * Routes to appropriate method based on k-mer size (small vs big).
	 * @param bases Sequence bases as byte array
	 * @param quals Quality scores (may be null)
	 */
	public final void hash(final byte[] bases, final byte[] quals){
		if(k<32){hashSmall(bases, quals);}
		else{hashBig(bases, quals);}
	}

	/**
	 * Hashes and tracks sequence using short k-mers (k < 32).
	 * Uses bit-shifting operations for efficient k-mer rolling.
	 * Applies quality-based probability filtering when quality scores provided.
	 * Uses canonical k-mers (maximum of forward and reverse complement).
	 * @param bases Sequence bases as byte array
	 * @param quals Quality scores for probability calculation (may be null)
	 */
	public final void hashSmall(final byte[] bases, final byte[] quals){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		final int lenmask=(-1)>>>k;

		long kmer=0, rkmer=0;

		if(minProb>0 && quals!=null){//Debranched loop
			assert(quals.length==bases.length) : quals.length+", "+bases.length;
			float prob=1;
			int len=0;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];

//				long x=AminoAcid.baseToNumber[b];
//				final long x2=x^3;

				final long x=((b>>1)&3)|((b&8)<<28);//Different encoding
				final long x2=x^2;

				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				//Update probability
				prob=prob*PROB_CORRECT[quals[i]];
				if(len>=k){prob=prob*PROB_CORRECT_INVERSE[quals[i-k]];}

				if(x>=0){
					len++;
				}else{
					len=0;
					kmer=rkmer=0;
					prob=1;
				}

				if(len>=k && prob>=minProb){add(Math.max(kmer, rkmer));}
			}
		}else{

			//Old, slower logic
//			for(int i=0; i<bases.length; i++){
//				byte b=bases[i];
//				long x=AminoAcid.baseToNumber[b];
//				long x2=AminoAcid.baseToComplementNumber[b];
//				kmer=((kmer<<2)|x)&mask;
//				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
//
//				if(x>=0){
//					len++;
//				}else{
//					len=0;
//					kmer=rkmer=0;
//				}
//				if(len>=k){
//					add(Math.max(kmer, rkmer));
//				}
//			}

			//Does not handle IUPAC
			int len=-1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
//				final long x=AminoAcid.baseToNumber[b];//Normal 2-bit encoding
//				final long x2=x^3;
				final long x=((b>>1)&3)|((b&8)<<28);//Alternate 2-bit encoding
				final long x2=x^2;

				// shift register, tracks valid bases in upper zeros
				len|=x;
				len>>>=1;

				// Let the kmers roll. Garbage flushes out naturally after k shifts!
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				// Countdown timer reached the threshold
				if(len<=lenmask){
					add(Math.max(kmer, rkmer));
				}
			}
		}
	}

	/**
	 * Hashes and tracks sequence using long k-mers (k >= 32).
	 * Uses Kmer objects for handling k-mers longer than 31 bases.
	 * Applies quality-based probability filtering when quality scores provided.
	 * @param bases Sequence bases as byte array
	 * @param quals Quality scores for probability calculation (may be null)
	 */
	public final void hashBig(final byte[] bases, final byte[] quals){
		final Kmer kmer=getLocalKmer();
		int len=0;
		float prob=1;

		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final long x=Dedupe.baseToNumber[b];
			kmer.addRightNumeric(x);
			if(minProb>0 && quals!=null){//Update probability
				prob=prob*PROB_CORRECT[quals[i]];
				if(len>=k){
					final byte oldq=quals[i-k];
					prob=prob*PROB_CORRECT_INVERSE[oldq];
				}
			}
			if(AminoAcid.isFullyDefined(b)){
				len++;
			}else{
				len=0;
				prob=1;
			}
			if(len>=k && prob>=minProb){
				add(kmer.xor());
			}
		}
	}

	/** Returns the compensation factor for the current bucket count. */
	public final float compensationFactorBuckets(){
		assert(Integer.bitCount(buckets)==1) : buckets;
		final int zeros=Integer.numberOfTrailingZeros(buckets);
		return compensationFactorLogBuckets(zeros);
	}

	/** Returns the compensation factor for the given log2(buckets) value. */
	public final float compensationFactorLogBuckets(final int logBuckets){
		final float[] array=compensationFactorLogBucketsArray();
		return (array!=null && logBuckets<array.length) ? array[logBuckets] : 1/(1+(1<<logBuckets));
	}

	/** Builds a sorted frequency list from bucket counts for histogram output. */
	public SuperLongList toFrequency(){
		final SuperLongList list=new SuperLongList(1000);
		if(getClass()==LogLog16.class){
			final char[] counts=counts16();
			for(final char x : counts){
				if(x>0){list.add(x);}
			}
		}else{
			final int[] counts=getCounts();
			for(final int x : counts){
				if(x>0){list.add(x);}
			}
		}
		list.sort();
		return list;
	}

	/**
	 * Prints a k-mer frequency histogram to file.
	 * Outputs depth-count pairs with optional supersampling adjustment.
	 * Handles both array-based and list-based frequency data.
	 * @param path File path for output
	 * @param overwrite Whether to overwrite existing file
	 * @param append Whether to append to existing file
	 * @param supersample Adjust counts for the effect of subsampling
	 * @param decimals Number of decimal places for supersampled counts
	 */
	/**
	 * Prints a k-mer frequency histogram to file.
	 * Outputs depth-count pairs with optional supersampling adjustment.
	 * @param path File path for output
	 * @param overwrite Whether to overwrite existing file
	 * @param append Whether to append to existing file
	 * @param supersample Adjust counts for the effect of subsampling
	 * @param decimals Number of decimal places for supersampled counts
	 */
	public void printKhist(final String path, final boolean overwrite, final boolean append, final boolean supersample, final int decimals){
		printKhist32(path, overwrite, append, supersample, decimals);
	}

	/** Prints a 32-bit k-mer frequency histogram to file. */
	public void printKhist32(final String path, final boolean overwrite, final boolean append, final boolean supersample, final int decimals){
		final SuperLongList sll=toFrequency();
		final ByteStreamWriter bsw=new ByteStreamWriter(path, overwrite, append, false);
		bsw.start();
		bsw.print("#Depth\tCount\n");
		final double mult=Math.max(1.0, (supersample ? cardinality()/(double)buckets : 1));
		final long[] array=sll.array();
		final LongList list=sll.list();

		for(int depth=0; depth<array.length; depth++){
			final long count=array[depth];
			if(count>0){
				bsw.print(depth).tab();
				if(supersample){
					if(decimals>0){
						bsw.print(count*mult, decimals).nl();
					}else{
						bsw.print(Math.max(1, Math.round(count*mult))).nl();
					}
				}else{
					bsw.print(count).nl();
				}
			}
		}
		int count=0;
		long prevDepth=-1;
		for(int i=0; i<list.size; i++){
			final long depth=list.get(i);
			if(depth!=prevDepth && count>0){
				assert(depth>prevDepth);
				bsw.print(prevDepth).tab();
				if(supersample){
					if(decimals>0){
						bsw.print(count*mult, decimals).nl();
					}else{
						bsw.print(Math.max(1, Math.round(count*mult))).nl();
					}
				}else{
					bsw.print(count).nl();
				}
				count=0;
			}else{
				count++;
			}
			prevDepth=depth;
		}
		if(count>0){
			bsw.print(prevDepth).tab();
			if(supersample){
				if(decimals>0){
					bsw.print(count*mult, decimals).nl();
				}else{
					bsw.print(Math.max(1, Math.round(count*mult))).nl();
				}
			}else{
				bsw.print(count).nl();
			}
		}
		bsw.poisonAndWait();
	}

//	public void printKhist16(String path, boolean overwrite, boolean append, boolean supersample, int decimals){
//		char[] array=counts16();
//		ByteStreamWriter bsw=new ByteStreamWriter(path, overwrite, append, false);
//		bsw.start();
//		bsw.print("#Depth\tCount\n");
//		final double mult=Math.max(1.0, (supersample ? cardinality()/(double)buckets : 1));
//		
//		for(int depth=0; depth<array.length; depth++){
//			long count=array[depth];
//			if(count>0){
//				bsw.print(depth).tab();
//				if(supersample){
//					if(decimals>0){
//						bsw.print(count*mult, decimals).nl();
//					}else{
//						bsw.print(Math.max(1, Math.round(count*mult))).nl();
//					}
//				}else{
//					bsw.print(count).nl();
//				}
//			}
//		}
//		bsw.poisonAndWait();
//	}

	/** Sum of all bucket counts. Returns 0 if counts are not tracked. */
	public final long countSum(){
		final int[] counts=getCounts();
		return counts==null ? 0 : simd.Vector.sum(counts);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Calculates cardinality estimate from this tracker.
	 * Implementation varies by subclass algorithm.
	 * @return Estimated number of unique k-mers observed
	 */
	public abstract long cardinality();

	/**
	 * Returns the counts array if present.
	 * Should be overridden for classes that track counts.
	 * @return Array of bucket counts, or null if not tracking counts
	 */
	public int[] getCounts(){
		return null;
	}

	/** Returns 16-bit counts array if present; null otherwise. */
	public char[] counts16(){return null;}
	public byte[] gcArray(){return null;}

	/**
	 * Merges another tracker into this one.
	 * Combines cardinality estimates from both trackers.
	 * @param log The tracker to add to this one
	 */
	public abstract void add(final CardinalityTracker log);

	/**
	 * Generates a 64-bit hashcode from a number and adds it to this tracker.
	 * <p>
	 * <b>DO NOT CALL DIRECTLY.</b> Use {@link #add(long)} instead, which calls
	 * this method and also increments the {@code added} counter needed for
	 * {@code clampToAdded} support. Calling hashAndStore directly will leave
	 * {@code added=0}, causing {@code cardinality()} to return 0 when clamping
	 * is enabled (the default).
	 * <p>
	 * This method exists as the subclass extension point — subclasses override
	 * this to implement their specific hashing and storage logic.
	 *
	 * @param number The value to hash and store
	 */
	public abstract void hashAndStore(final long number);

	/**
	 * Returns array of compensation factors indexed by log2(buckets).
	 * Designed to compensate for overestimate with small numbers of buckets.
	 * May be deprecated in favor of harmonic mean handling multiple trials.
	 * @return Array of compensation factors, or null if not implemented
	 */
	public abstract float[] compensationFactorLogBucketsArray();

	/**
	 * Terminal correction for the plain Mean estimator.
	 * At high cardinality the uncorrected Mean converges to some class-specific
	 * ratio of the true cardinality; this method returns that ratio so it can
	 * be divided out before CF-table lookup. Default 1.0 leaves Mean unchanged.
	 * Override in subclasses with the empirically measured asymptotic bias
	 * (average CF_Mean across the top octave of a preliminary CF table).
	 * When the override matches reality, regenerated CF tables converge to 1.0.
	 */
	public float terminalMeanCF(){return 1f;}

	/**
	 * Terminal correction for the history-corrected Mean+H estimator.
	 * Only meaningful for classes with history bits. Default 1.0 leaves Mean+H
	 * unchanged. Override in UDLL6 and other history classes with the
	 * empirically measured asymptotic bias of Mean+H.
	 */
	public float terminalMeanPlusCF(){return 1f;}

	/** HC weight for LDLC blend: LDLC = w*HC + (1-w)*DLC.
	 *  Default 0.50 = equal blend of HC and DLC.
	 *  Override in subclasses with calibrated values. */
	public double ldlcHcWeight(){return 0.50;}

	/** LDLC weight in HLDLC blend: HLDLC = w*LDLC + (1-w)*Hybrid+2.
	 *  Default 0.5; override per class based on empirical sweep. */
	public float hldlcWeight(){return OVERRIDE_HLDLC_WEIGHT>=0 ? OVERRIDE_HLDLC_WEIGHT : 0.5f;}

	/*--------------------------------------------------------------*/
	/*----------------       Drivable Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Default: throws UnsupportedOperationException. Override in calibratable subclasses. */
	@Override
	public double[] rawEstimates(){throw new UnsupportedOperationException(getClass().getSimpleName()+" does not support rawEstimates()");}

	/** Default: throws UnsupportedOperationException. Override in calibratable subclasses. */
	@Override
	public int filledBuckets(){throw new UnsupportedOperationException(getClass().getSimpleName()+" does not support filledBuckets()");}

	/** Default: throws UnsupportedOperationException. Override in calibratable subclasses. */
	@Override
	public double occupancy(){throw new UnsupportedOperationException(getClass().getSimpleName()+" does not support occupancy()");}

	@Override
	public long getLastCardinality(){return lastCardinality;}

	@Override
	public void setLastCardinality(final long val){lastCardinality=val;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** K-mer length for sequence hashing. */
	public final int k;

	/** Minimum probability threshold for k-mer correctness filtering. */
	public final float minProb;

	/** Number of buckets for tracking; must be a power of 2 for efficiency. */
	public final int buckets;
	/** Returns the actual bucket count used by this tracker. */
	public int actualBuckets(){return buckets;}

	/** Log2 of buckets, i.e. number of bits in bucketMask. */
	public final int bucketBits;

	/** Bit mask for extracting bucket index from hashcode; equals buckets-1. */
	final int bucketMask;

	/** Thread-local storage for Kmer objects used in long k-mer hashing. */
	private final ThreadLocal<Kmer> localKmer=new ThreadLocal<Kmer>();

	/** Value XORed with hash inputs to vary hash function behavior. */
	protected final long hashXor;

	/**
	 * Gets a thread-local Kmer object for long k-mer mode processing.
	 * Creates new Kmer if none exists for current thread.
	 * Clears the Kmer before returning for reuse.
	 * @return Thread-local Kmer object ready for use
	 */
	protected Kmer getLocalKmer(){
		Kmer kmer=localKmer.get();
		if(kmer==null){
			synchronized(localKmer){
				kmer=localKmer.get();
				if(kmer==null){
					kmer=new Kmer(k);
					localKmer.set(kmer);
				}
			}
		}
		kmer.clearFast();
		return kmer;
	}

	/** Number of hashes added to this tracker. */
	long added=0;
	/** 64-bit micro cardinality sketch for very-low-cardinality estimation. */
	long microIndex=0;
	/** Absolute NLZ histogram: nlzCounts[k] = number of buckets with absoluteNlz==k.
	 * Lazy-allocated on first summarize() call; cleared and reused on subsequent calls.
	 * Used by DynamicLC (DLC) estimator in CardinalityStats. Null for non-DLL/DDL classes. */
	protected int[] nlzCounts;
	/** Sorted dif-value buffer for median, MWA, and trimmed-mean estimators.
	 * Lazy-allocated on first summarize() when USE_SORTBUF is true; null otherwise.
	 * Gated by USE_SORTBUF so HotSpot eliminates the dead code when disabled. */
	protected LongList sortBuf;
	/** Cached cardinality estimate; -1 means stale. Shared by all calibratable subclasses. */
	public long lastCardinality=-1;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Array converting quality scores to probability of base correctness. */
	public static final float[] PROB_CORRECT=Arrays.copyOf(align2.QualityTools.PROB_CORRECT, 128);
	/** Inverse probability array: 1/PROB_CORRECT for sliding-window probability updates. */
	public static final float[] PROB_CORRECT_INVERSE=Arrays.copyOf(align2.QualityTools.PROB_CORRECT_INVERSE, 128);

	// --- Global configuration ---
	// These fields are set once by CardinalityParser.initializeAll() on the main thread
	// before any worker threads are created. Thread.start() provides happens-before
	// visibility to worker threads. Do not mutate after threads are running.
	//TODO: Replace with a final immutable config container published before thread creation

	public static boolean trackCounts=false;
	public static boolean clampToAdded=true;
	public static long lastCardinalityStatic=-1;
	public static long defaultSeed=0;
	public static boolean USE_MEAN=true;
	public static boolean USE_MEDIAN=false;
	public static boolean USE_MWA=false;
	public static boolean USE_HMEAN=false;
	public static boolean USE_HMEANM=false;
	public static boolean USE_GMEAN=false;
	public static boolean USE_HLL=false;
	public static boolean USE_HYBRID=false;
	public static final boolean LAZY_ALLOCATE=false;
	/** When true, allocate and fill sortBuf for median/MWA/mean99 estimators.
	 * When false (default), those estimators return fallback values and no LongList is allocated. */
	public static final boolean USE_SORTBUF=false;
	/** When true, use cardinality-indexed CF table in addition to occupancy-indexed table.
	 * Only applies when CorrectionFactor.USE_CORRECTION is also true; cardinality CF
	 * cannot be active without the occupancy CF. Default true. */
	public static boolean USE_CARD_CF=true;
	/** When true, use history bits to improve LC at low cardinality.
	 * Each non-empty bucket contributes 1 + popcount(history bits) to a
	 * set-bit total, providing a tighter lower bound on cardinality. */
	public static boolean USE_HISTORY_FOR_LC=false;
	/** HC blend weight for LDLC estimator.
	 * Set from ldlcHcWeight() instance method in initializeAll(). */
	public static double LDLC_HC_WEIGHT=0.50;
	/** Command-line override for ldlcHcWeight(). Negative = use class default. */
	public static double OVERRIDE_LDLC_HC_WEIGHT=-1;
	/** Command-line override for hldlcWeight(). Negative = use class default. */
	public static float OVERRIDE_HLDLC_WEIGHT=-1f;
	/** When true AND SBS table is loaded, hybrid estimators use sbs()
	 * instead of lcMin as the low-cardinality component in the blend zone. */
	public static boolean USE_SBS_IN_HYBRID=true;
	public static final int MICRO_CUTOFF_BITS=56;//Higher is less accurate, max is 64

	/*--------------------------------------------------------------*/
	/*----------------    Reflective Factory         ----------------*/
	/*--------------------------------------------------------------*/

	private static CardinalityTracker makeByReflection(String cls){
		try{return (CardinalityTracker)Class.forName("cardinality."+cls).getDeclaredConstructor().newInstance();}
		catch(Exception e){throw new RuntimeException(cls, e);}
	}

	private static CardinalityTracker makeByReflection(String cls, parse.Parser p){
		try{return (CardinalityTracker)Class.forName("cardinality."+cls).getDeclaredConstructor(parse.Parser.class).newInstance(p);}
		catch(Exception e){throw new RuntimeException(cls, e);}
	}

	private static CardinalityTracker makeByReflection(String cls, int buckets, int k, long seed, float minProb){
		try{return (CardinalityTracker)Class.forName("cardinality."+cls).getDeclaredConstructor(int.class, int.class, long.class, float.class).newInstance(buckets, k, seed, minProb);}
		catch(Exception e){throw new RuntimeException(cls, e);}
	}

}
