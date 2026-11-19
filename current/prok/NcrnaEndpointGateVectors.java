package prok;

import java.io.File;
import java.util.Locale;

import fileIO.ReadWrite;
import map.LongHashSet;
import map.ObjectIntMap;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Makes second-stage move/hold rows from frozen first-stage predictions.
 * Fold mode permits only explicitly held internal assemblies. Outer mode reuses
 * the already saved baseline scores and never evaluates the first-stage net.
 * Both modes keep the independent-end trial assay: invalid proposed spans are
 * reported, not silently removed or relabeled. Actual caller joint holds remain
 * a separate integration requirement. Nothing here fits or selects a gate.
 * Adapted from NcrnaPositionalFeatureVariants and EvaluateNcrnaPositionalNet.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaEndpointGateVectors{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointGateVectors.class, false);
		try{
			new NcrnaEndpointGateVectors(pp.args).process();
			pp.outstream.println("NCRNA_ENDPOINT_GATE_VECTORS_COMPLETE_NO_GATE_FIT");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaEndpointGateVectors(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("rows")){rows=b;}
			else if(a.equals("net")){network=b;}
			else if(a.equals("scores")){saved=b;}
			else if(a.equals("contexts")){contexts=b;}
			else if(a.equals("ref")){reference=b;}
			else if(a.equals("original")){original=b;}
			else if(a.equals("graph")){graph=b;}
			else if(a.equals("roster")){roster=b;}
			else if(a.equals("loci")){loci=b;}
			else if(a.equals("exclude")){exclude=b;}
			else if(a.equals("mapping")){mapping=b;}
			else if(a.equals("fold")){fold=Integer.parseInt(b);}
			else if(a.equals("samples")){samples=Integer.parseInt(b);}
			else if(a.equals("expected")){expected=Long.parseLong(b);}
			else if(a.equals("minlen")){minimumLength=Integer.parseInt(b);}
			else if(a.equals("end")){
				require(b.equals("5prime") || b.equals("3prime"), "end= must identify one endpoint");
				end=b.equals("5prime") ? 0 : 1;
			}else{require(parser.parse(arg, a, b), "Unknown gate-vector argument: "+arg);}
		}
		matrix=parser.in1;
		output=parser.out1;
		require(matrix!=null && output!=null && rows!=null && network!=null && contexts!=null && reference!=null
			&& original!=null && graph!=null && roster!=null && loci!=null && exclude!=null && mapping!=null,
			"Bind in= rows= net= contexts= ref= original= graph= roster= loci= exclude= mapping= out=fresh_directory");
		require(end>=0 && samples>0 && expected>0 && minimumLength>0 && fold>=-1 && fold<5,
			"Require end= samples= expected=accepted_producer_count minlen= and fold=-1(outer) or0..4");
		require((fold==-1)==(saved!=null), "Outer rows must reuse scores=; fold rows require one new first-stage inference");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1 && !new File(output).exists(),
			"Gate rows require fresh outputs and complete input populations, without sampling");
		require(Tools.testInputFiles(false, true, matrix, rows, network, contexts, reference, original, graph, roster, loci, exclude, mapping),
			"Every gate-row population, model and context binding must be readable");
		require(saved==null || Tools.testInputFiles(false, true, saved), "Saved first-stage score stream must be readable");
	}

	public void process(){process(null);}

	/** The optional borrowed network exists for finite fixtures; CLI always loads its bound file. */
	public void process(final CellNet fixtureNet){
		require(fixtureNet==null || saved==null, "Saved outer scores must never be replaced by network inference");
		population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			NcrnaSupport.loadLibrary(original).names, roster, loci, exclude);
		library=NcrnaSupport.loadLibrary(reference);
		loadMapping();
		loadContexts();
		modelGC=new float[library.bases.length];
		for(int i=0; i<modelGC.length; i++){
			modelGC[i]=NcrnaPositionalFeatureVariants.gc(library.bases[i], 0, library.bases[i].length-1);
		}
		final CellNet net=saved==null ? fixtureNet==null ? CellNetParser.load(network, false) : fixtureNet : null;
		require(net==null || net.numInputs()==69 && net.numOutputs()==22, "Fold nets must retain the frozen69/22recipe");
		require(new File(output).mkdir(), "Cannot create fresh gate-vector output directory");
		final boolean bgzf=ReadWrite.USE_BGZF;
		final boolean force=ReadWrite.FORCE_BGZIP;
		final boolean explicit=ReadWrite.SET_ZIP_THREADS;
		final int threads=ReadWrite.MAX_ZIP_THREADS();
		ReadWrite.USE_BGZF=ReadWrite.FORCE_BGZIP=true;
		ReadWrite.setZipThreads(1);
		try{stream(net);}
		finally{
			ReadWrite.USE_BGZF=bgzf;
			ReadWrite.FORCE_BGZIP=force;
			ReadWrite.setZipThreads(threads);
			ReadWrite.SET_ZIP_THREADS=explicit;
		}
	}

	/** Mapping must cover the entire original axis and agree on every exact assembly. */
	private void loadMapping(){
		assigned=new int[population.ids.length];
		java.util.Arrays.fill(assigned, -2);
		final ObjectIntMap<String> assemblies=new ObjectIntMap<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(mapping)){
			RrnaResourceIO.header(input.reader, "locus_id\tsource\tassembly\trole\tfold\tpartial", mapping);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==6, "Fold mapping needs all six fields");
				final int gene=population.index(p.parseString(0));
				require(gene>=0 && assigned[gene]==-2 && p.termEquals(population.sourceName(gene), 1)
					&& p.termEquals(population.assemblyName(gene), 2), "Fold mapping must identify each original exactly once");
				final int role=p.parseInt(4);
				require(role>=-1 && role<5 && p.termEquals(role==-1 ? "OUTER_TEST" : "INTERNAL_FOLD", 3)
					&& p.termEquals(population.partial[gene] ? "1" : "0", 5)
					&& population.excluded[gene]==(role==-1), "Outer exclusions and original partial dispositions cannot change");
				final String assembly=population.assemblyName(gene);
				final int old=assemblies.get(assembly);
				require(old<0 || old==role+1, "An exact assembly cannot occupy two fold roles");
				assemblies.put(assembly, role+1);
				assigned[gene]=role;
			}
		}
		for(final int role:assigned){require(role!=-2, "Fold mapping must conserve every original locus");}
	}

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
				require(bases[gene]==null && read.bases!=null, "Each original needs one nonempty context");
				bases[gene]=read.bases;
				contigGC[gene]=context.gc;
				available[gene]=context.hasFeatures();
			});
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
		}
		for(final byte[] sequence:bases){require(sequence!=null, "An original context is missing");}
	}

	/** Streaming joins preserve identities; all predictions are saved before gating. */
	private void stream(final CellNet net){
		final LineParser1 p=new LineParser1('\t');
		final LineParser1 meta=new LineParser1('\t');
		final LineParser1 scoresRow=new LineParser1('\t');
		final LongHashSet seen=new LongHashSet();
		final float[] features=new float[69];
		final float[] scores=new float[22];
		final float[] gate=new float[32];
		final ByteBuilder text=new ByteBuilder();
		long total=0;
		long proposed=0;
		long positives=0;
		long crossed=0;
		long shortSpan=0;
		long invalidPositive=0;
		int last=-1;
		try(Input data=new Input(matrix); Input identities=new Input(rows);
				Input replay=saved==null ? null : new Input(saved);
				Output vectors=new Output(path("vectors.tsv.gz")); Output audit=new Output(path("rows.tsv"));
				Output raw=saved==null ? new Output(path("first_stage.tsv.gz")) : null){
			RrnaResourceIO.header(data.reader, "#dims\t69\t22\t0", matrix);
			RrnaResourceIO.header(identities.reader, NcrnaPositionalFeatureVariants.ROW_HEADER, rows);
			if(replay!=null){readScoreHeader(replay);}
			vectors.writer.println("#dims\t32\t1\t0");
			audit.writer.println(NcrnaPositionalFeatureVariants.ROW_HEADER+
				"\tgate_row1\tproposal_class\ttrue_minus_raw_bp\tstrictly_better\tcrossed\tbelow_minimum\toor_head_above_proposal");
			if(raw!=null){
				raw.writer.print(new ByteBuilder("#NCRNA_GATE_FIRST_STAGE_V1\n#coding\tFP32_RAW_SIGNED_INT\n")
					.append("#matrix\t").append(matrix).nl().append("#rows\t").append(rows).nl()
					.append("#network\t").append(network).nl().append("#end\t").append(end==0 ? "5prime" : "3prime").nl());
				text.clear().append(NcrnaPositionalFeatureVariants.ROW_HEADER);
				for(int i=0; i<22; i++){text.append("\tscore_bits_").append(i);}
				raw.writer.print(text.nl());
			}
			for(byte[] line=data.reader.nextLine(); line!=null; line=data.reader.nextLine()){
				final byte[] identity=identities.reader.nextLine();
				require(identity!=null, "Gate sidecar ended before its matrix");
				p.set(line);
				meta.set(identity);
				require(p.terms()==91 && meta.terms()==11, "Gate inputs require the complete69/22matrix and eleven-field sidecar");
				for(int i=0; i<11; i++){require(meta.length(i)>0, "Gate sidecar fields cannot be empty");}
				for(int i=0; i<69; i++){
					features[i]=number(p, i);
					require(Float.isFinite(features[i]), "Nonfinite first-stage input");
				}
				final int gene=last>=0 && meta.termEquals(population.ids[last], 0) ? last : population.index(meta.parseString(0));
				require(gene>=0 && assigned[gene]==fold && available[gene] && !population.partial[gene]
					&& meta.termEquals(population.sourceName(gene), 1) && meta.termEquals(population.assemblyName(gene), 2)
					&& meta.termEquals(fold==-1 ? "evaluation" : "FOLD_HELD", 3), "Gate row must belong to its declared held assembly role");
				last=gene;
				final int sample=meta.parseInt(6);
				final int model=meta.parseInt(4);
				require(sample>=0 && sample<samples && seen.add(((long)gene<<32)|sample), "Duplicate or invalid held trial");
				require(model>=0 && model<library.names.length && meta.termEquals(library.names[model], 5), "Selected model differs from its frozen catalog");
				final int start=meta.parseInt(7);
				final int stop=meta.parseInt(8);
				final int truthStart=population.leftFlank[gene];
				final int truthStop=bases[gene].length-population.rightFlank[gene]-1;
				require(start>=0 && stop>=start && stop<bases[gene].length && meta.parseInt(9)==truthStart
					&& meta.parseInt(10)==truthStop, "Raw and true coordinates must bind the original context");
				final long correction=end==0 ? (long)truthStart-start : (long)truthStop-stop;
				final int label=correction<-10 || correction>10 ? 21 : (int)correction+10;
				for(int i=0; i<22; i++){require(p.termEquals(i==label ? "1" : "0", 69+i), "First-stage target disagrees with exact endpoint correction");}
				final int refLength=library.bases[model].length;
				require(features[63]==(stop-start+1)/(float)refLength && features[64]==contigGC[gene]
					&& features[65]>=0 && features[65]<=1, "Raw length, GC or identity disagrees with original feature bindings");
				if(replay==null){
					net.applyInput(features);
					net.feedForward();
					for(int i=0; i<22; i++){scores[i]=net.getOutput(i);}
				}else{
					readSavedRow(replay, scoresRow, identity, label, correction, scores);
				}
				final int choice=proposal(scores);
				if(raw!=null){
					text.clear().append(identity);
					for(final float score:scores){text.tab().append(Float.floatToRawIntBits(score));}
					raw.writer.print(text.nl());
				}
				total++;
				boolean better=false;
				boolean cross=false;
				boolean tooShort=false;
				if(choice>=0){
					final int shift=choice-10;
					final long length=(long)stop-start+1+(end==0 ? -shift : shift);
					cross=length<=0;
					tooShort=length<minimumLength;
					better=Math.abs(correction-shift)<Math.abs(correction);
					fill(scores, features, NcrnaPositionalFeatureVariants.gc(bases[gene], start, stop),
						modelGC[model], length, refLength, gate);
					text.clear();
					for(final float value:gate){text.append(value, 18).tab();}
					vectors.writer.print(text.append(better ? 1 : 0).nl());
					proposed++;
					if(better){positives++;}
					if(cross){crossed++;}
					if(tooShort){shortSpan++;}
					if(tooShort && better){invalidPositive++;}
				}
				audit.writer.print(text.clear().append(identity).tab().append(choice<0 ? 0 : proposed).tab().append(choice)
					.tab().append(correction).tab().append(better ? 1 : 0).tab().append(cross ? 1 : 0)
					.tab().append(tooShort ? 1 : 0).tab().append(choice>=0 && scores[21]>scores[choice] ? 1 : 0).nl());
			}
			require(total==expected && identities.reader.nextLine()==null, "Gate rows must conserve the accepted producer count: expected="+expected+" actual="+total);
			if(replay!=null){
				final byte[] terminal=replay.reader.nextLine();
				require(terminal!=null, "Saved scores require their terminal completion marker");
				scoresRow.set(terminal);
				require(scoresRow.terms()==2 && scoresRow.termEquals("#COMPLETE", 0)
					&& scoresRow.parseLong(1)==total && replay.reader.nextLine()==null, "Saved scores must close at the same complete trial count");
			}else{raw.writer.print(text.clear().append("#COMPLETE\t").append(total).nl());}
		}
		try(Output summary=new Output(path("summary.tsv"))){
			summary.writer.println("fold\tend\ttrials\tproposals\tpositives\tholds\tcrossed\tbelow_minimum\tinvalid_positives_retained\tinputs\tfirst_stage_inference\tfeature_order");
			summary.writer.print(text.clear().append(fold).tab().append(end==0 ? "5prime" : "3prime").tab().append(total)
				.tab().append(proposed).tab().append(positives).tab().append(total-proposed).tab().append(crossed).tab().append(shortSpan)
				.tab().append(invalidPositive).append("\t32\t").append(saved==null ? "ONCE" : "SAVED_ONLY").tab().append(FEATURE_ORDER).nl());
		}
	}

	/** A tie or the current position winning bypasses the gate. Head21 is a feature, not a veto. */
	public static int proposal(final float[] scores){return NcrnaEndpointGateFeatures.proposal(scores);}

	/** Length features compare original/proposed spans; GC, identity and fuzz remain raw-span features. */
	public static void fill(final float[] scores, final float[] raw, final float spanGC, final float refGC,
			final long proposedLength, final int referenceLength, final float[] result){
		NcrnaEndpointGateFeatures.fill(scores, raw, spanGC, refGC, proposedLength, referenceLength, result);
	}

	private void readScoreHeader(final Input input){
		RrnaResourceIO.header(input.reader, "#NCRNA_FULL_ERROR_SCORES_V1", saved);
		RrnaResourceIO.header(input.reader, "#coding\tFP32_RAW_SIGNED_INT", saved);
		RrnaResourceIO.header(input.reader, "#matrix\t"+matrix, saved);
		RrnaResourceIO.header(input.reader, "#rows\t"+rows, saved);
		RrnaResourceIO.header(input.reader, "#network\t"+network, saved);
		RrnaResourceIO.header(input.reader, "#end\t"+(end==0 ? "5prime" : "3prime"), saved);
		byte[] line=input.reader.nextLine();
		final LineParser1 p=new LineParser1('\t');
		require(line!=null, "Saved scores lack their column header");
		p.set(line);
		if(p.termEquals("#featurevariant", 0)){
			require(p.terms()==2 && p.termEquals("baseline", 1), "Gate comparison reuses only the frozen baseline first stage");
			line=input.reader.nextLine();
		}
		final ByteBuilder header=new ByteBuilder(NcrnaPositionalFeatureVariants.ROW_HEADER)
			.append("\ttrue_class\ttrue_minus_raw_bp\tmargin0_class\tmargin01_class\targmax22_class");
		for(int i=0; i<22; i++){header.append("\tscore_bits_").append(i);}
		require(line!=null && java.util.Arrays.equals(line, header.toBytes()), "Saved score column order must match the qualified evaluator");
	}

	private static void readSavedRow(final Input input, final LineParser1 p, final byte[] identity,
			final int label, final long correction, final float[] scores){
		final byte[] line=input.reader.nextLine();
		require(line!=null, "Saved first-stage scores ended before the paired input");
		p.set(line);
		require(p.terms()==38, "Saved score rows require the qualified38-field schema");
		p.setBounds(10);
		require(p.b()==identity.length, "Saved score identity width differs from its exact sidecar");
		for(int i=0; i<identity.length; i++){require(line[i]==identity[i], "Saved score and sidecar identities disagree");}
		require(p.parseInt(11)==label && p.parseLong(12)==correction, "Saved score truth must equal the paired coordinates");
		for(int i=0; i<22; i++){
			require(p.length(16+i)>0, "Raw score bits cannot be empty");
			scores[i]=Float.intBitsToFloat(p.parseInt(16+i));
			require(Float.isFinite(scores[i]), "Saved first-stage scores must be finite");
		}
		// The legacy endpoint-window argmax requires odd width. This finite
		// 22-head schema uses the native unrestricted first-maximum helper.
		require(p.parseInt(13)==NcrnaEndpointConfidence.selectFull(scores, 10, 0)
			&& p.parseInt(14)==NcrnaEndpointConfidence.selectFull(scores, 10, .1f)
			&& p.parseInt(15)==Tools.maxIndex(scores), "Saved policy choices must agree with preserved scores");
	}
	private static float number(final LineParser1 p, final int term){
		require(p.length(term)>0, "Numeric feature fields cannot be empty");
		p.setBounds(term);
		for(int i=p.a(); i<p.b(); i++){
			if(p.line()[i]=='e' || p.line()[i]=='E'){return (float)Parse.parseDoubleSlow(p.line(), p.a(), p.b());}
		}
		return p.parseFloat(term);
	}
	private String path(final String name){return new File(output, name).getPath();}
	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}
	public static final String FEATURE_ORDER=NcrnaEndpointGateFeatures.FEATURE_ORDER;
	private String matrix, rows, network, saved, contexts, reference, original, graph, roster, loci, exclude, mapping, output;
	private int fold=-2, samples=-1, minimumLength=-1, end=-1;
	private long expected=-1;
	private NcrnaLearningPopulation population;
	private NcrnaSupport.Library library;
	private int[] assigned;
	private byte[][] bases;
	private float[] contigGC, modelGC;
	private boolean[] available;
}
