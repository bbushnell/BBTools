package var2;

import stream.Read;

/**
 * Regression tests for merging FASTA scaffold names with SAM/VCF short names.
 */
public class ScafMapNameTest {

	public static void main(String[] args){
		testVcfMergesWithReference();
		testSamMergesWithReference();
		testReferenceMergesWithPriorSam();
		testReferenceRejectsPriorSamLengthCollision();
		testSameLengthShortNameReuse();
		testDistinctFullNamesRemainDistinct();
		testExactShortNamePrecedence();
		System.out.println("ScafMap name tests passed.");
	}

	private static void testVcfMergesWithReference(){
		ScafMap map=new ScafMap();
		Scaffold ref=map.addScaffold(read("chr1 descriptive text", "ACGT"));
		Scaffold vcf=map.addFromVcf("##contig=<ID=chr1,length=4>".getBytes());
		check(ref==vcf, "VCF short name created a duplicate scaffold");
		check(map.size()==1, "VCF merge changed scaffold count");
		check(map.getScaffold("chr1").bases!=null, "VCF lookup lost reference bases");
	}

	private static void testSamMergesWithReference(){
		ScafMap map=new ScafMap();
		Scaffold ref=map.addScaffold(read("chr1 descriptive text", "ACGT"));
		Scaffold sam=map.add("@SQ\tSN:chr1\tLN:4".getBytes());
		check(ref==sam, "SAM short name created a duplicate scaffold");
		check(map.size()==1, "SAM merge changed scaffold count");
	}

	private static void testReferenceMergesWithPriorSam(){
		ScafMap map=new ScafMap();
		Scaffold sam=map.add("@SQ\tSN:chr1\tLN:4".getBytes());
		Scaffold ref=map.addScaffold(read("chr1 descriptive text", "ACGT"));
		check(ref==sam, "reference created a duplicate of a prior SAM scaffold");
		check(map.size()==1, "reverse-order SAM/reference merge changed scaffold count");
		check(map.getScaffold("chr1")==sam, "short lookup failed after reference merge");
		check(map.getScaffold("chr1 descriptive text")==sam, "full lookup failed after reference merge");
		check(sam.bases!=null, "reference bases were not attached to the SAM scaffold");
	}

	private static void testReferenceRejectsPriorSamLengthCollision(){
		ScafMap map=new ScafMap();
		map.add("@SQ\tSN:chr1\tLN:5".getBytes());
		try{
			map.addScaffold(read("chr1 descriptive text", "ACGT"));
			throw new AssertionError("Expected reverse-order length collision");
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains("ambiguous"), "wrong length-collision exception");
		}
		check(map.size()==1, "length collision mutated scaffold count");
	}

	private static void testSameLengthShortNameReuse(){
		ScafMap map=new ScafMap();
		Scaffold sam=map.add("@SQ\tSN:dup\tLN:4".getBytes());
		Scaffold first=map.addScaffold(read("dup first", "AAAA"));
		Scaffold second=map.addScaffold(read("dup second", "AAAA"));
		check(first==sam && second==sam, "equal-length short-name scaffolds were not reused");
		check(map.size()==1, "equal-length short-name reuse changed scaffold count");
	}

	private static void testDistinctFullNamesRemainDistinct(){
		ScafMap map=new ScafMap();
		Scaffold first=map.addScaffold(read("dup first", "AAAA"));
		Scaffold second=map.addScaffold(read("dup second", "CCCC"));
		check(first!=second, "distinct full names were collapsed");
		check(map.size()==2, "distinct full names changed scaffold count");
		check(map.getScaffold("dup first")==first, "first exact lookup failed");
		check(map.getScaffold("dup second")==second, "second exact lookup failed");
		try{
			map.addFromVcf("##contig=<ID=dup,length=4>".getBytes());
			throw new AssertionError("Expected ambiguous VCF short-name failure");
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains("Ambiguous"), "wrong ambiguity exception");
		}
	}

	private static void testExactShortNamePrecedence(){
		ScafMap map=new ScafMap();
		Scaffold exact=map.addScaffold(read("dup", "AAAA"));
		Scaffold longer=map.addScaffold(read("dup second", "CCCC"));
		Scaffold vcf=map.addFromVcf("##contig=<ID=dup,length=4>".getBytes());
		check(vcf==exact, "exact short name did not take precedence");
		check(map.getScaffold("dup second")==longer, "long exact name was lost");
		check(map.size()==2, "exact short-name lookup changed scaffold count");
	}

	private static Read read(String name, String bases){
		return new Read(bases.getBytes(), null, name, 0);
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
