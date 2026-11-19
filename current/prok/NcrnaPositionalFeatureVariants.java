package prok;

import java.io.File;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.Outputs;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Projects three independent feature changes from the same retained trials.
 * Copies unchanged numeric fields verbatim, never realigns, redraws, recounts
 * endpoint tables or fits a network. Exact original-locus, assembly, context,
 * model and target joins guard against plausible but mismatched sidecars.
 * Each output keeps its input sidecar by reference, not by a bulk duplicate.
 * Adapted from the native positional-vector and saved-evaluation byte readers.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaPositionalFeatureVariants{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPositionalFeatureVariants.class, false);
		try{
			new NcrnaPositionalFeatureVariants(pp.args).process();
			pp.outstream.println("NCRNA_POSITIONAL_FEATURE_VARIANTS_COMPLETE");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaPositionalFeatureVariants(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("contexts")){contexts=b;}
			else if(a.equals("ref")){reference=b;}
			else if(a.equals("original")){original=b;}
			else if(a.equals("graph")){graph=b;}
			else if(a.equals("roster")){roster=b;}
			else if(a.equals("loci")){loci=b;}
			else if(a.equals("exclude")){exclusions=b;}
			else if(a.equals("samples")){samples=Integer.parseInt(b);}
			else if(a.equals("expected")){
				final String[] counts=b.split(",", -1);
				require(counts.length==6, "expected= binds train5,train3,monitor5,monitor3,evaluation5,evaluation3 counts");
				expected=new long[6];
				for(int i=0; i<6; i++){
					expected[i]=Long.parseLong(counts[i]);
					require(expected[i]>0, "Every retained matrix must contain trials");
				}
			}else{require(parser.parse(arg, a, b), "Unknown option: "+arg);}
		}
		input=parser.in1;
		output=parser.out1;
		require(input!=null && output!=null && contexts!=null && reference!=null && original!=null
			&& graph!=null && roster!=null && loci!=null && exclusions!=null && expected!=null && samples>0,
			"Require in=matrix_directory out=fresh_directory contexts= ref= original= graph= roster= loci= exclude= samples= expected=six_counts");
		require(parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1, "Project complete retained inputs to fresh outputs without sampling");
		require(new File(input).isDirectory() && !new File(output).exists(), "Input directory must exist and output directory must be fresh");
		require(Tools.testInputFiles(false, true, contexts, reference, original, graph, roster, loci, exclusions),
			"All exact population/context/model bindings must be readable");
		for(int set=0; set<3; set++){
			for(int end=0; end<2; end++){
				final String stem=stem(set, end);
				require(Tools.testInputFiles(false, true, path(input, stem+".tsv"), path(input, stem+".rows.tsv")),
					"Each of the six matrices needs its matching row sidecar");
			}
		}
	}

	public void process(){
		population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			NcrnaSupport.loadLibrary(original).names, roster, loci, exclusions);
		library=NcrnaSupport.loadLibrary(reference);
		modelGC=new float[library.bases.length];
		minimum=Integer.MAX_VALUE;
		maximum=0;
		for(int model=0; model<library.bases.length; model++){
			final int length=library.bases[model].length;
			require(length>0, "Reference length normalization needs nonempty models");
			minimum=Math.min(minimum, length);
			maximum=Math.max(maximum, length);
			modelGC[model]=gc(library.bases[model], 0, length-1);
		}
		validateAssemblies();
		loadContexts();
		require(new File(output).mkdir(), "Cannot create fresh projection output directory");
		for(final String variant:VARIANTS){require(new File(output, variant).mkdir(), "Cannot create fresh variant directory");}
		final boolean bgzf=ReadWrite.USE_BGZF;
		final boolean force=ReadWrite.FORCE_BGZIP;
		final boolean zipExplicit=ReadWrite.SET_ZIP_THREADS;
		final int zipThreads=ReadWrite.MAX_ZIP_THREADS();
		ReadWrite.USE_BGZF=ReadWrite.FORCE_BGZIP=true;
		ReadWrite.setZipThreads(1);
		try(Output summary=new Output(path(output, "projection.tsv"))){
			summary.writer.println("variant\tsplit\tend\trows\tinputs\toutputs\tref_min\tref_max\tquery_scaled_above_one\tidentity_scaled_below_zero\tmatrix\tunchanged_sidecar\tfeature_order");
			for(int set=0; set<3; set++){
				for(int end=0; end<2; end++){project(set, end, summary.writer);}
			}
		}finally{
			ReadWrite.USE_BGZF=bgzf;
			ReadWrite.FORCE_BGZIP=force;
			ReadWrite.setZipThreads(zipThreads);
			ReadWrite.SET_ZIP_THREADS=zipExplicit;
		}
	}

	/** Exclusions act on complete genomes even though the saved file names loci. */
	private void validateAssemblies(){
		final ObjectIntMap<String> membership=new ObjectIntMap<String>(String.class);
		for(int gene=0; gene<population.ids.length; gene++){
			if(population.partial[gene]){continue;}
			final String assembly=population.assemblyName(gene);
			require(!assembly.equals(".") && !assembly.equals("NA"), "Every retained trial requires exact assembly identity");
			final int split=population.excluded[gene] ? 1 : 0;
			final int previous=membership.get(assembly);
			require(previous<0 || previous==split, "One assembly cannot teach and evaluate: "+assembly);
			membership.put(assembly, split);
		}
	}

	/** Contexts are small locus sequences, retained once; no genome FASTA is opened. */
	private void loadContexts(){
		bases=new byte[population.ids.length][];
		contigGC=new float[bases.length];
		available=new boolean[bases.length];
		final boolean trim=Shared.TRIM_READ_DESCRIPTION;
		final boolean upper=Read.TO_UPPER_CASE;
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
		final NcrnaLearningContext context=new NcrnaLearningContext();
		try{
			RrnaResourceIO.stream(contexts, read->{
				context.set(read, population);
				final int gene=context.gene;
				require(bases[gene]==null && read.bases!=null, "Contexts must contain one nonempty record per original locus");
				bases[gene]=read.bases;
				contigGC[gene]=context.gc;
				available[gene]=context.hasFeatures();
			});
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
		}
		for(int gene=0; gene<bases.length; gene++){require(bases[gene]!=null, "Missing original context: "+population.ids[gene]);}
	}

	/** One streaming pass supplies all three single-variable alternatives. */
	private void project(final int set, final int end, final ByteStreamWriter summary){
		final String stem=stem(set, end);
		final String sidecar=path(input, stem+".rows.tsv");
		final LineParser1 values=new LineParser1('\t');
		final LineParser1 metadata=new LineParser1('\t');
		final float[] parsed=new float[69];
		final ByteBuilder[] buffers={new ByteBuilder(), new ByteBuilder(), new ByteBuilder()};
		final int[] seen=new int[Math.multiplyExact(population.ids.length, samples)];
		long count=0;
		long above=0;
		long below=0;
		int last=-1;
		try(Input matrix=new Input(path(input, stem+".tsv")); Input rows=new Input(sidecar); Outputs outputs=new Outputs(3)){
			RrnaResourceIO.header(matrix.reader, "#dims\t69\t22\t0", stem);
			RrnaResourceIO.header(rows.reader, ROW_HEADER, sidecar);
			final ByteStreamWriter[] writers=new ByteStreamWriter[3];
			for(int variant=0; variant<3; variant++){
				// This single producer writes sequentially through print(); ordered
				// ByteStreamWriter mode instead requires numbered add() jobs.
				writers[variant]=outputs.open(FileFormat.testOutput(destination(variant, stem), FileFormat.TEXT, null, true, false, false, false));
				writers[variant].print(new ByteBuilder("#dims\t").append(WIDTHS[variant]).append("\t22\t0\n"));
			}
			for(byte[] line=matrix.reader.nextLine(); line!=null; line=matrix.reader.nextLine()){
				final byte[] row=rows.reader.nextLine();
				require(row!=null, "Sidecar ended before matrix: "+stem);
				values.set(line);
				metadata.set(row);
				require(values.terms()==91 && metadata.terms()==11, "Full-error matrix/sidecar field counts changed: "+stem);
				for(int i=0; i<11; i++){require(metadata.length(i)>0, "Sidecar fields must be explicit, never implicit zero values");}
				for(int i=0; i<91; i++){require(values.length(i)>0, "Empty numeric fields cannot define retained trials");}
				for(int i=0; i<69; i++){
					parsed[i]=numericFloat(values, i);
					require(Float.isFinite(parsed[i]), "Nonfinite retained input: "+stem);
				}
				final int gene=last>=0 && metadata.termEquals(population.ids[last], 0) ? last : population.index(metadata.parseString(0));
				require(gene>=0 && available[gene] && !population.partial[gene], "Sidecar requires an available complete original context");
				last=gene;
				require(metadata.termEquals(population.sourceName(gene), 1) && metadata.termEquals(population.assemblyName(gene), 2)
					&& metadata.termEquals(SETS[set], 3) && population.excluded[gene]==(set==2), "Row identity/split differs from its original population");
				final int model=metadata.parseInt(4);
				final int sample=metadata.parseInt(6);
				require(model>=0 && model<library.names.length && metadata.termEquals(library.names[model], 5), "Selected model differs from the ordered reference witness");
				require(sample>=0 && sample<samples && seen[gene*samples+sample]++==0, "Duplicate or invalid locus/sample in retained matrix");
				final int start=metadata.parseInt(7);
				final int stop=metadata.parseInt(8);
				final int trueStart=population.leftFlank[gene];
				final int trueStop=bases[gene].length-population.rightFlank[gene]-1;
				require(metadata.parseInt(9)==trueStart && metadata.parseInt(10)==trueStop
					&& start>=0 && stop>=start && stop<bases[gene].length, "Retained raw/truth span disagrees with context geometry");
				final long correction=end==0 ? (long)trueStart-start : (long)trueStop-stop;
				final int label=correction<-10 || correction>10 ? 21 : (int)correction+10;
				for(int i=0; i<22; i++){
					require(values.termEquals(i==label ? "1" : "0", 69+i), "One-hot label disagrees with the exact signed endpoint correction");
				}
				final int length=library.bases[model].length;
				require(parsed[63]==(stop-start+1)/(float)length && parsed[64]==contigGC[gene]
					&& parsed[65]>=0 && parsed[65]<=1, "Original length/contig GC/identity inputs disagree with their bindings");
				final float refLength=referenceLength(length, minimum, maximum);
				final float queryLength=(stop-start+1)/(2f*length);
				final float spanGC=gc(bases[gene], start, stop);
				final float identity=2f*parsed[65]-1f;
				if(queryLength>1){above++;}
				if(identity<0){below++;}
				projectRow(values, refLength, queryLength, spanGC, modelGC[model], identity, buffers);
				for(int variant=0; variant<3; variant++){
					writers[variant].print(buffers[variant]);
					buffers[variant].clear();
				}
				count++;
			}
			require(rows.reader.nextLine()==null && count==expected[2*set+end], "Projection must conserve every bound matrix/sidecar row: "+stem+" rows="+count);
		}
		for(int variant=0; variant<3; variant++){
			summary.print(new ByteBuilder().append(VARIANTS[variant]).tab().append(SETS[set]).tab().append(ENDS[end])
				.tab().append(count).tab().append(WIDTHS[variant]).append("\t22\t").append(minimum).tab().append(maximum)
				.tab().append(above).tab().append(below).tab().append(destination(variant, stem)).tab().append(sidecar)
				.tab().append(ORDERS[variant]).nl());
		}
	}

	/** Same retained numeric contract as EvaluateNcrnaPositionalNet: ordinary
	 * decimals stay on the native fast path; valid e/E tokens use its slow parser.
	 */
	private static float numericFloat(final LineParser1 parser, final int term){
		parser.setBounds(term);
		final byte[] line=parser.line();
		for(int i=parser.a(); i<parser.b(); i++){
			if(line[i]=='e' || line[i]=='E'){return (float)Parse.parseDoubleSlow(line, parser.a(), parser.b());}
		}
		return parser.parseFloat(term);
	}

	/** Finite fixture entrypoint; unchanged fields and targets retain exact bytes. */
	public static void projectRow(final LineParser1 row, final float refLength, final float queryLength,
			final float spanGC, final float referenceGC, final float identity, final ByteBuilder[] output){
		require(row!=null && row.terms()==91 && output!=null && output.length==3,
			"Feature projection requires one complete 69/22 row and three owned output buffers");
		copy(row, output[0], 0, 63);
		output[0].append(refLength, 18).tab().append(queryLength, 18).tab();
		copy(row, output[0], 64, 91);
		copy(row, output[1], 0, 65);
		output[1].append(spanGC, 18).tab().append(referenceGC, 18).tab();
		copy(row, output[1], 65, 91);
		copy(row, output[2], 0, 65);
		output[2].append(identity, 18).tab();
		copy(row, output[2], 66, 91);
		for(final ByteBuilder buffer:output){buffer.nl();}
	}

	/** Copies [first,last) fields, including the separating tab when not final. */
	private static void copy(final LineParser1 row, final ByteBuilder out, final int first, final int last){
		assert(first>=0 && last>first && last<=row.terms()) : "Projection copies only complete bounded original fields";
		row.setBounds(first);
		final int from=row.a();
		row.setBounds(last-1);
		out.append(row.line(), from, row.b()-from);
		if(last<row.terms()){out.tab();}
	}

	/** GC/(A+C+G+T), matching Read.gc's recognized-base denominator; all-N is unavailable. */
	public static float gc(final byte[] sequence, final int start, final int stop){
		return NcrnaEndpointGateFeatures.gc(sequence, start, stop);
	}

	public static float referenceLength(final int length, final int minimum, final int maximum){
		require(minimum>0 && minimum<=length && length<=maximum, "Reference length must lie inside its frozen library range");
		return minimum==maximum ? 0 : (length-minimum)/(float)(maximum-minimum);
	}
	private String destination(final int variant, final String stem){return path(path(output, VARIANTS[variant]), stem+".tsv.gz");}
	private static String stem(final int set, final int end){return "r10."+SETS[set]+'.'+ENDS[end];}
	private static String path(final String directory, final String file){return new File(directory, file).getPath();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	public static final String ROW_HEADER="locus_id\tsource\tassembly\tsplit\tmodel_index\tmodel\tsample\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0";
	private static final String[] SETS={"train", "monitor", "evaluation"};
	private static final String[] ENDS={"5prime", "3prime"};
	private static final String[] VARIANTS={"length", "gc", "identity"};
	private static final int[] WIDTHS={70, 71, 69};
	private static final String PREFIX="bestScores,leftScores,rightScores,";
	private static final String FUZZ=",hbmCoverageMinus,hbmCoverageAt,hbmCoveragePlus";
	private static final String[] ORDERS={PREFIX+"referenceLength,queryLength,contigGC,candidateIdentity"+FUZZ,
		PREFIX+"lengthRatio,contigGC,spanGC,referenceGC,candidateIdentity"+FUZZ,
		PREFIX+"lengthRatio,contigGC,rescaledIdentity"+FUZZ};
	private String input, output, contexts, reference, original, graph, roster, loci, exclusions;
	private long[] expected;
	private int samples;
	private NcrnaLearningPopulation population;
	private NcrnaSupport.Library library;
	private byte[][] bases;
	private float[] contigGC, modelGC;
	private boolean[] available;
	private int minimum, maximum;
}
