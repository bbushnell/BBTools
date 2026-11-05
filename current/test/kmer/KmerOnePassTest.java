package test.kmer;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;

/**
 * Finite one-pass count and quality contracts using six identical three-base reads.
 * A depth-two prefilter intentionally withholds the first two accepted occurrences;
 * the all-high main table should retain four. Mixed cases independently gate filter
 * support and main counting, and check input counts before depth rejection.
 * Run each case/assertion mode in a fresh JVM and a separate retained input directory.
 * @author Sucrose
 */
public final class KmerOnePassTest {

	/** Selects one independent count/quality case and its retained input directory. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("dir=")){dir=Paths.get(arg.substring(4));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(mode==null || dir==null){throw new IllegalArgumentException("Expected case=<case name> and dir=<new input directory>");}
		final boolean onePass, lowQuality, filterQuality, mainQuality, checkInputCount;
		final int expectedStored;
		String qualityPattern=null;
		int expectedInput=6;
		if(mode.equals("count-onepass")){
			onePass=true; lowQuality=false; filterQuality=true; mainQuality=false; checkInputCount=true; expectedStored=4;
		}else if(mode.equals("count-unfiltered")){
			onePass=false; lowQuality=false; filterQuality=false; mainQuality=true; checkInputCount=true; expectedStored=6;
		}else if(mode.equals("quality-prefilter")){
			onePass=true; lowQuality=true; filterQuality=true; mainQuality=false; checkInputCount=false; expectedStored=0;
		}else if(mode.equals("quality-main")){
			onePass=true; lowQuality=true; filterQuality=false; mainQuality=true; checkInputCount=false; expectedStored=0;
		}else if(mode.equals("quality-both")){
			onePass=true; lowQuality=true; filterQuality=true; mainQuality=true; checkInputCount=false; expectedStored=0;
		}else if(mode.equals("quality-disabled")){
			onePass=true; lowQuality=true; filterQuality=false; mainQuality=false; checkInputCount=false; expectedStored=4;
		}else if(mode.equals("quality-high")){
			onePass=true; lowQuality=false; filterQuality=true; mainQuality=false; checkInputCount=false; expectedStored=4;
		}else if(mode.equals("quality-unfiltered")){
			onePass=false; lowQuality=true; filterQuality=false; mainQuality=true; checkInputCount=false; expectedStored=0;
		}else if(mode.equals("mixed-low-first-none")){
			onePass=true; lowQuality=false; filterQuality=false; mainQuality=false; checkInputCount=true; expectedStored=4;
			qualityPattern="LLLHHH";
		}else if(mode.equals("mixed-low-first-main")){
			onePass=true; lowQuality=false; filterQuality=false; mainQuality=true; checkInputCount=true; expectedStored=3;
			qualityPattern="LLLHHH"; expectedInput=3;
		}else if(mode.equals("mixed-low-first-prefilter")){
			onePass=true; lowQuality=false; filterQuality=true; mainQuality=false; checkInputCount=true; expectedStored=1;
			qualityPattern="LLLHHH";
		}else if(mode.equals("mixed-low-first-both")){
			onePass=true; lowQuality=false; filterQuality=true; mainQuality=true; checkInputCount=true; expectedStored=1;
			qualityPattern="LLLHHH"; expectedInput=3;
		}else if(mode.equals("mixed-high-first-none")){
			onePass=true; lowQuality=false; filterQuality=false; mainQuality=false; checkInputCount=true; expectedStored=4;
			qualityPattern="HHHLLL";
		}else if(mode.equals("mixed-high-first-main")){
			onePass=true; lowQuality=false; filterQuality=false; mainQuality=true; checkInputCount=true; expectedStored=1;
			qualityPattern="HHHLLL"; expectedInput=3;
		}else if(mode.equals("mixed-high-first-prefilter")){
			onePass=true; lowQuality=false; filterQuality=true; mainQuality=false; checkInputCount=true; expectedStored=4;
			qualityPattern="HHHLLL";
		}else if(mode.equals("mixed-high-first-both")){
			onePass=true; lowQuality=false; filterQuality=true; mainQuality=true; checkInputCount=true; expectedStored=1;
			qualityPattern="HHHLLL"; expectedInput=3;
		}else if(mode.equals("mixed-boundary")){
			onePass=true; lowQuality=false; filterQuality=true; mainQuality=false; checkInputCount=true; expectedStored=0;
			qualityPattern="HHLLLL";// Two qualifying observations cannot pass the depth-two filter.
		}else{throw new IllegalArgumentException("Unknown case: "+mode);}

		Files.createDirectories(dir);
		final Path input=dir.resolve("reads.fq");
		final StringBuilder fastq=new StringBuilder();
		for(int i=0; i<6; i++){
			final boolean low=qualityPattern==null ? lowQuality : qualityPattern.charAt(i)=='L';
			fastq.append("@read").append(i).append("\nACG\n+\n").append(low ? "###\n" : "III\n");
		}
		Files.write(input, fastq.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		AbstractKmerTableSet.MASK_CORE=false;
		AbstractKmerTableSet.MASK_MIDDLE=false;
		final KmerTableSet tables=new KmerTableSet(new String[] {
			"in="+input, "k=3", "ways=1", "threads=1", "prealloc=f",
			"prefilter="+(onePass ? "2" : "0"), "onepass="+onePass, "prepasses=1", "filtermem=1m",
			"minprob=0.9", "minprobprefilter="+filterQuality, "minprobmain="+mainQuality,
			"rcomp=t", "qin=33", "interleaved=f", "minavgquality=0", "qtrim=f", "ecco=f", "merge=f"
		}, 12);
		check(tables.onePass==onePass && tables.prefilter==onePass, "Unexpected requested counting mode");
		tables.processInput();
		check(!tables.errorState, "Input loading reported an error");
		check(tables.prefilter==onePass, "Prefilter was unexpectedly disabled for lack of memory");
		check(tables.readsIn==6 && tables.basesIn==18, "Expected six untrimmed three-base reads");
		check(tables.lowqReads==0 && tables.readsTrimmed==0, "Whole-read filtering must not mask k-mer probability filtering");
		check(tables.kmersLoaded==(expectedStored>0 ? 1 : 0), "Wrong number of distinct retained k-mers: "+tables.kmersLoaded);
		for(int key=0; key<64; key++){
			final int expected=key==27 ? expectedStored : 0;// ACG's canonical reverse complement CGT encodes27.
			check(Math.max(0, tables.getCount(key))==expected,
				"Wrong retained count for key "+key+": "+tables.getCount(key)+" versus "+expected);
		}
		if(checkInputCount){
			check(tables.kmersIn==expectedInput, "Input count must apply main quality filtering before depth rejection: "+
				tables.kmersIn+" versus "+expectedInput);
		}
		System.out.println("PASS case="+mode+" reads=6 stored="+expectedStored+" checkedInputCount="+checkInputCount);
	}

	/** Keeps count/quality oracles active in both Java assertion modes. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
