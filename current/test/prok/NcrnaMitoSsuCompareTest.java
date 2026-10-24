package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import idaligner.AlignmentStats;
import prok.NcrnaMitoSsuCompare;
import prok.NcrnaMitoSsuCompare.Result;
import structures.ByteBuilder;

/** Explicit identity/orientation oracles plus full CLI multiplicity and parallel parity.
 * @author Brian Bushnell, Raiden */
public final class NcrnaMitoSsuCompareTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="), "Require outdir=NEW_DIRECTORY");
		final String root=args[0].substring(7);require(new File(root).mkdir(), "Fixture output directory must be new");
		final String a="ACGTGACTGACCTGAACTGA", c="TTGCAACCTTAGGCTAACGC", fragment=a.substring(4, 14), longA="TTTT"+a+"CCCC";
		final Result exact=compare(a, a), shortGene=compare(fragment, a), longGene=compare(longA, a);
		require(exact.identity==1f && exact.traceIdentity==1f && exact.matches==a.length() && exact.candidateQuery, "Identical equal-length input is a complete exact candidate-query alignment");
		require(shortGene.identity==1f && shortGene.paired==fragment.length() && shortGene.start==4 && shortGene.stop==13 && shortGene.candidateQuery,
			"Exact internal fragment must retain model-side coordinates and partial model coverage");
		require(longGene.identity==1f && longGene.paired==a.length() && longGene.start==4 && longGene.stop==23 && !longGene.candidateQuery,
			"When model is shorter, raw endpoints belong to the candidate and paired fraction excludes its flanks");
		final Result neutral=compare("ACGTNACGT", "ACGTAACGT");
		require(neutral.ns==1 && neutral.matches==8 && Math.abs(neutral.identity-8.5f/9)<0.00001f && Math.abs(neutral.traceIdentity-8f/9)<0.00001f,
			"Quantum returned N half credit and trace zero credit are distinct documented conventions");
		write(root+"/ref16.fa", ">universal\n"+a+"\n>unused\n"+c+"\n");
		write(root+"/ref18.fa", ">euk18S_universal\n"+c+"\n");
		final String input=">one metadata\n"+a+"\n>two\n"+a+"\n>three\n"+c+"\n>fragment\n"+fragment+"\n>long\n"+longA+"\n>six\n"+a.toLowerCase(java.util.Locale.ROOT)+"\n";
		write(root+"/in.fa", input);run(root, "serial", root+"/in.fa", root+"/ref16.fa", 6, 1);run(root, "parallel", root+"/in.fa", root+"/ref16.fa", 6, 3);
		for(String suffix:new String[]{".pairs.tsv.gz", ".members.tsv.gz", ".unique.fa.gz", ".summary.tsv"}){
			require(read(root+"/serial"+suffix).equals(read(root+"/parallel"+suffix)), "Thread count must not change dedup IDs, alignments, member order or summaries: "+suffix);
		}
		final String pairs=read(root+"/serial.pairs.tsv.gz"), members=read(root+"/serial.members.tsv.gz"), summary=read(root+"/serial.summary.tsv");
		require(lines(pairs)==9 && lines(members)==7 && lines(summary)==3, "Four unique sequences need eight alignments and all six member mappings");
		require(members.contains("one\tseq_0\t20\n") && members.contains("two\tseq_0\t20\n") && members.contains("six\tseq_0\t20\n"),
			"Exact sequence caching must preserve all duplicate copies, normalizing input case only");
		require(summary.contains("copy\t6\t5\t0\t1\t") && summary.contains("unique\t4\t3\t0\t1\t"), "Copy-weighted and unique-weighted preference denominators must differ honestly");
		final ByteBuilder many=new ByteBuilder();
		for(int i=0; i<40; i++){
			many.append('>').append("many_").append(i).nl().append(a);
			for(int shift=0; shift<8; shift+=2){many.append("ACGT".charAt((i>>>shift)&3));}many.nl();
		}
		write(root+"/many.fa", many.toString());run(root, "many_serial", root+"/many.fa", root+"/ref16.fa", 40, 1);run(root, "many_parallel", root+"/many.fa", root+"/ref16.fa", 40, 3);
		for(String suffix:new String[]{".pairs.tsv.gz", ".members.tsv.gz", ".unique.fa.gz", ".summary.tsv"}){
			require(read(root+"/many_serial"+suffix).equals(read(root+"/many_parallel"+suffix)), "At least three16-item chunks must preserve serial/parallel output parity: "+suffix);
		}
		require(lines(read(root+"/many_parallel.pairs.tsv.gz"))==81, "Forty distinct sequences must produce eighty alignments across multiple work chunks");
		write(root+"/duplicate.fa", input+">one\n"+a+"\n");reject(root, "duplicate", root+"/duplicate.fa", root+"/ref16.fa", 7, "Duplicate source locus");
		reject(root, "count", root+"/in.fa", root+"/ref16.fa", 7, "Input locus count differs");
		write(root+"/missing_ref.fa", ">other\n"+a+"\n");reject(root, "missing_ref", root+"/in.fa", root+"/missing_ref.fa", 6, "Missing exact universal model");
		write(root+"/duplicate_ref.fa", ">universal\n"+a+"\n>universal\n"+c+"\n");reject(root, "duplicate_ref", root+"/in.fa", root+"/duplicate_ref.fa", 6, "Duplicate universal model");
		System.out.println("MITO_SSU_COMPARE_TEST_PASS loci=6 unique=4 alignments=8 parallel_unique=40 negatives=4 parallel_exact=true fragments=true query_swap=true neutral_N=true");
	}
	private static Result compare(String a, String b){return NcrnaMitoSsuCompare.compare(a.getBytes(StandardCharsets.US_ASCII), b.getBytes(StandardCharsets.US_ASCII), new AlignmentStats(true));}
	private static void run(String root, String name, String input, String ref, int expected, int threads){
		NcrnaMitoSsuCompare.main(new String[]{"in="+input, "ref16="+ref, "ref18="+root+"/ref18.fa", "outprefix="+root+"/"+name, "expected="+expected, "t="+threads, "zl=9"});
	}
	private static void reject(String root, String name, String input, String ref, int expected, String message){
		boolean failed=false;try{run(root, name, input, ref, expected, 1);}catch(IllegalArgumentException e){require(e.getMessage().contains(message), "Unexpected failure: "+e);failed=true;}
		require(failed && !new File(root+"/"+name+".summary.tsv").exists(), "Invalid input cannot publish a success summary: "+name);
	}
	private static int lines(String s){int n=0;for(int i=0; i<s.length(); i++){if(s.charAt(i)=='\n'){n++;}}return n;}
	private static String read(String path){
		final ByteFile f=ByteFile.makeByteFile(path, true);final ByteBuilder b=new ByteBuilder();
		try{for(byte[] line; (line=f.nextLine())!=null;){b.append(line).nl();}}finally{require(!f.close(), "Fixture read failed: "+path);}
		return new String(b.array, 0, b.length(), StandardCharsets.UTF_8);
	}
	private static void write(String path, String value){
		final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path, FileFormat.TEXT, null, true, false, false, false));w.start();w.print(value.getBytes(StandardCharsets.UTF_8));require(!w.poisonAndWait(), "Fixture write failed");
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
}
