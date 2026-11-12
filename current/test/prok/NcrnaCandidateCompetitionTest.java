package test.prok;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.Data;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference;
import prok.NcrnaCallerProcess;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;

/** Full-caller path-competition regression on the unchanged Microvirgula
 * aerodenitrificans DSM15089 assembly GCF_000620105.1. Bulk input stays external.
 * Adapted from template.A_SampleBasic2; native launchers own genomic I/O/workers.
 * A FIXED adversarial configuration disallows same-strand overlap, boosts TPP's
 * path score and excludes shortened CDS starts. Only the shortlist cap changes.
 * Explicit target/locus arguments can bind another independently identified fixture.
 * Ordinary controls must prove the named RNA/CDS acceptance flip BEFORE reuse
 * executes. Fresh reuse-off/on outputs must equal both complete ordinary GFFs
 * and native reference grades. This configuration is never a production choice.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCandidateCompetitionTest{

	public static void main(final String[] args){
		final NcrnaCandidateCompetitionTest test=new NcrnaCandidateCompetitionTest(args);
		try{test.process();}
		finally{Shared.closeStream(test.outstream);}
	}

	public NcrnaCandidateCompetitionTest(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCandidateCompetitionTest.class, false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(java.util.Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("reference")){reference=b;}
			else if(a.equals("familymap")){familymap=b;}
			else if(a.equals("target")){target=b;}
			else if(a.equals("gradefamily")){gradeFamily=b;}
			else if(a.equals("contig")){
				contig=b;
				coordinates|=1;
			}else if(a.equals("strand")){
				strand=b;
				coordinates|=2;
			}else if(a.equals("rnastart")){
				rnaStart=Integer.parseInt(b);
				coordinates|=4;
			}else if(a.equals("rnastop")){
				rnaStop=Integer.parseInt(b);
				coordinates|=8;
			}else if(a.equals("cdsstart")){
				cdsStart=Integer.parseInt(b);
				coordinates|=16;
			}else if(a.equals("cdsstop")){
				cdsStop=Integer.parseInt(b);
				coordinates|=32;
			}else if(a.equals("mincdslen")){minCdsLength=Integer.parseInt(b);}
			else if(a.equals("baselinecap")){baselineCap=Integer.parseInt(b);}
			else if(a.equals("floor")){floor=Float.parseFloat(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		check(parser.in1!=null && parser.out1!=null && reference!=null && familymap!=null,
			"Supply in=original_genome reference=registered_manifest familymap=rules out=fresh_directory");
		check("tpp".equals(target) || "sam1".equals(target), "This competition fixture currently supports TPP or SAM-I");
		check(target.equals("tpp") || coordinates==63 && "SAM_I".equals(gradeFamily),
			"A SAM-I fixture needs an explicit grading family and all independently identified locus coordinates");
		check(contig!=null && !contig.isEmpty() && ("+".equals(strand) || "-".equals(strand)) && rnaStart>0 && rnaStop>=rnaStart
			&& cdsStart>0 && cdsStop>=cdsStart && minCdsLength>0 && baselineCap>1
			&& Float.isFinite(floor) && floor>=0, "Fixture coordinates and policy must be valid");
		genome=new File(parser.in1).getAbsolutePath();
		root=new File(parser.out1).getAbsoluteFile();
		check(Tools.testInputFiles(false, true, genome, reference, familymap), "Fixture inputs must be readable");
		for(final String path:new String[]{genome, reference, familymap, root.getPath()}){
			check(path.indexOf('\n')<0 && path.indexOf('\r')<0 && path.indexOf('\t')<0,
				"Native configuration and literal manifest paths must occupy one field");
		}
		check(!root.exists() && root.mkdir(), "Use a fresh fixture directory under an existing parent");
	}

	private void process(){
		String inputId=null;
		for(final GradeNcrnaReference.Binding binding:GradeNcrnaReference.loadManifest(reference)){
			if(new File(binding.fasta).getAbsolutePath().equals(genome)){
				check(inputId==null, "The original fixture genome must have exactly one reference binding");
				inputId=binding.id;
			}
		}
		check(inputId!=null, "The fixture genome must be present in the supplied registered reference");
		final String[] common={"in="+genome, "taxonomy=f", "callcds=t", "calltrna=t", "call16s=t",
			"call23s=t", "call5s=t", "call18s=f", "ncrna=t", "tmrna=t", "sixs=t", "fmn=t",
			"sam1=t", "tpp=t", "tppcdsstop=t", "glycine=t", "cobalamin=t", "euk=f",
			"ncrnafamily="+target, "ncrnascorea=100000", "maxoverlapss=0", "minlen="+minCdsLength,
			"ordered=t", "t=1", "ow=f", "append=f"};
		for(int p=0; p<2; p++){
			write(file("ordinary_"+p+".config"), common,
				"ncrnaindextopn="+(p==0 ? baselineCap : 1), "ncrnaadaptfloor="+floor,
				"out="+file("ordinary_"+p+".gff"));
			NcrnaCallerProcess.runInstalled(file("ordinary_"+p+".config"), 4096);
		}
		final int[] broad=features(file("ordinary_0.gff")), strict=features(file("ordinary_1.gff"));
		check(broad[0]==1 && broad[1]==0 && strict[0]==0 && strict[1]==1,
			"Ordinary controls must flip "+target+"/CDS acceptance; broad="+Arrays.toString(broad)+", strict="+Arrays.toString(strict));
		write(file("ordinary_flip.tsv"), new String[]{"policy\t"+gradeFamily+"_at_"+rnaStart+"_"+rnaStop+"\tCDS_at_"+cdsStart+"_"+cdsStop,
			"broad\t1\t0", "strict\t0\t1"});

		final long[][] ordinary=new long[2][];
		for(int p=0; p<2; p++){
			write(file("queries_"+p+".tsv"), new String[]{"input_id\tfasta\tgff",
				inputId+"\t"+genome+"\t"+file("ordinary_"+p+".gff")});
			write(file("grade_"+p+".config"), new String[]{"in="+file("queries_"+p+".tsv"),
				"reference="+reference, "familymap="+familymap, "restrictreference=t", "out="+file("grade_"+p)});
			run("dev/gradencrnareference.sh", file("grade_"+p+".config"), 1024);
			ordinary[p]=metrics(file("grade_"+p+".summary.tsv"), gradeFamily);
		}
		write(file("caller.config"), common);
		write(file("policies.tsv"), new String[]{"policy\ttopn\tfloor", "broad\t"+baselineCap+"\t"+floor, "strict\t1\t"+floor});
		for(final String mode:new String[]{"recomputed", "reused"}){
			write(file(mode+".config"), new String[]{"caller="+file("caller.config"), "reference="+reference,
				"familymap="+familymap, "target="+target, "policies="+file("policies.tsv"), "out="+file(mode),
				"savegff=t", "reusecandidates="+mode.equals("reused")});
			run("dev/sweepncrnashortlist.sh", file(mode+".config"), 4096);
			for(int p=0; p<2; p++){
				equalGff(file("ordinary_"+p+".gff"), file(mode+".input0.policy"+p+".gff"));
				check(Arrays.equals(ordinary[p], metrics(file(mode+".summary.tsv"), p==0 ? "broad" : "strict")),
					"Full-caller and native saved-GFF grades must agree for "+mode+" policy "+p);
			}
		}
		final String result="ORDINARY_CALLER_"+gradeFamily+"_CDS_FLIP_AND_REUSE_EQUIVALENCE_PASS";
		write(file("result.txt"), new String[]{result});
		outstream.println(result);
	}

	/** Count exact named genomic features; this is a regression assertion, not a grader. */
	private int[] features(final File path){
		final int[] result=new int[2];
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(path.getPath())){
			for(byte[] line; (line=next(input))!=null;){
				parser.set(line);
				check(parser.terms()==9, "Expected nine-column native GFF");
				if(!parser.termEquals(contig, 0) || !parser.termEquals(strand, 6)){continue;}
				final int start=parser.parseInt(3), stop=parser.parseInt(4);
				if(parser.termEquals("RNA", 2) && start==rnaStart && stop==rnaStop && Tools.indexOf(line, ",family:"+prok.NcrnaFamilyRegistry.displayName(target), 0)>=0){result[0]++;}
				if(parser.termEquals("CDS", 2) && start==cdsStart && stop==cdsStop){result[1]++;}
			}
		}
		return result;
	}

	public static void equalGff(final File expected, final File actual){
		try(Input a=new Input(expected.getPath()); Input b=new Input(actual.getPath())){
			long row=0;
			for(;;){
				final byte[] left=next(a), right=next(b);
				check(Arrays.equals(left, right), "Native GFF disagreement at record "+(++row)+": "+expected+" versus "+actual);
				if(left==null){return;}
			}
		}
	}

	public static long[] metrics(final File path, final String key){
		return metrics(path, key, null);
	}

	/** Optional second key selects one policy within a per-input sweep table. */
	public static long[] metrics(final File path, final String key, final String policy){
		final String[] names={"TP", "FP", "FN", "scored_references", "masked_references", "masked_calls",
			"exact_5prime", "exact_3prime", "exact_both", "duplicate_queries"};
		final int[] columns=new int[names.length];
		Arrays.fill(columns, -1);
		final LineParser1 parser=new LineParser1('\t');
		long[] result=null;
		try(Input input=new Input(path.getPath())){
			final byte[] header=next(input);
			check(header!=null, "Missing native grade/sweep header: "+path);
			parser.set(header);
			for(int i=0; i<names.length; i++){
				for(int c=0; c<parser.terms(); c++){
					if(parser.termEquals(names[i], c) || names[i].equals("duplicate_queries") && parser.termEquals("duplicates", c)){columns[i]=c;}
				}
				check(columns[i]>=0, "Missing metric "+names[i]+" in "+path);
			}
			for(byte[] line; (line=next(input))!=null;){
				parser.set(line);
				if(!parser.termEquals(key, 0) || policy!=null && !parser.termEquals(policy, 1)){continue;}
				check(result==null, "Metric row must be unique: "+key+" in "+path);
				result=new long[names.length];
				for(int i=0; i<columns.length; i++){result[i]=parser.parseLong(columns[i]);}
			}
		}
		check(result!=null, "Missing metric row "+key+" in "+path);
		return result;
	}

	private static byte[] next(final Input input){
		for(byte[] line; (line=input.reader.nextLine())!=null;){
			if(line.length>0 && line[0]!='#'){return line;}
		}
		return null;
	}

	/** Startup configuration strings may contain Unicode paths; encode explicitly. */
	public static void write(final File file, final String[] lines, final String... tail){
		try(Output output=new Output(file.getPath())){
			for(final String line:lines){output.writer.println(line.getBytes(StandardCharsets.UTF_8));}
			for(final String line:tail){output.writer.println(line.getBytes(StandardCharsets.UTF_8));}
		}
	}

	public static void run(final String script, final File config, final long heap){
		final String path=Data.findPath("?"+script, false);
		check(path!=null, "Missing native launcher "+script);
		NcrnaCallerProcess.run(new File(path), config, heap);
	}

	private File file(final String name){return new File(root, name);}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
	private final PrintStream outstream;
	private final String genome;
	private final File root;
	private String reference, familymap;
	private String target="tpp", gradeFamily="TPP", contig="NZ_JHVK01000002.1", strand="-";
	private int rnaStart=100514, rnaStop=100612, cdsStart=99648, cdsStop=100520;
	private int minCdsLength=870, baselineCap=361, coordinates;
	private float floor=5;
}
