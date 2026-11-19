package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import idaligner.EndClippedAligner;
import map.IntObjectMap;
import map.LongObjectMap;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;
import structures.LongList;
import template.ThreadWaiter;

/** Raw all-model measurements on one accepted capture shard, never genomes.
 * Adapted from A_SampleBasic2/A_SampleByteFileMT and RankSavedNcrnaWindows:
 * native byte metadata prepasses, Streamer FASTA input, worker-owned kernel
 * scratch, coarse claims and a joined owner reduction. Each retained window is
 * aligned once/model/stage before its observations are shared across genes.
 * Eligibility curves are NOT final caller recall, verifier replay parity or a
 * selected setting. Primary, physical-edge and their union remain explicit.
 * @author Brian Bushnell, Raiden
 */
public final class MeasureNcrnaRejectedWindows{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, MeasureNcrnaRejectedWindows.class, false);
		final Timer timer=new Timer();
		try{
			new MeasureNcrnaRejectedWindows(pp.args).process();
			timer.stop();
			pp.outstream.println("NCRNA_REJECTED_SAVED_WINDOW_MEASUREMENT_COMPLETE time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private MeasureNcrnaRejectedWindows(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("library")){label=b;}
			else if(a.equals("population")){population=b;}
			else if(a.equals("grade")){grade=b;}
			else if(a.equals("capture")){capture=b;}
			else if(a.equals("captureid")){captureId=b;}
			else if(a.equals("consensus")){consensus=b;}
			else if(a.equals("expectedmodels")){expectedModels=Parse.parseIntKMG(b);}
			else if(a.equals("expectedclass2")){expectedClass2=Parse.parseIntKMG(b);}
			else if(a.equals("expectedgenes")){expectedGenes=Parse.parseIntKMG(b);}
			else if(a.equals("expectedincomplete")){expectedIncomplete=Parse.parseIntKMG(b);}
			else if(a.equals("expectedtruncated")){expectedTruncated=Parse.parseIntKMG(b);}
			else if(a.equals("expectedshortboth")){expectedShortBoth=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown saved-attempt parameter: "+arg);}
		}
		misses=parser.in1;
		prefix=parser.out1;
		require(label!=null && label.matches("[A-Za-z0-9_]+") && captureId!=null && captureId.matches("[A-Za-z0-9_]+")
			&& population!=null && grade!=null && capture!=null && consensus!=null && misses!=null && prefix!=null
			&& expectedModels>0 && expectedClass2>0 && expectedGenes>0 && expectedIncomplete>=0 && expectedTruncated>=0 && expectedShortBoth>=0
			&& !parser.overwrite && !parser.append
			&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Supply exact library/population/grade/capture/captureid/consensus/in=misses/out=fresh paths and expected counts; no sampling");
		RrnaResourceIO.fresh(prefix, ".best.tsv.gz", ".curves.tsv", ".histograms.tsv", ".summary.tsv");
		capturePrefix=(captureId+":").getBytes(StandardCharsets.US_ASCII);
		workers=Shared.threads();
		require(workers>0 && workers<=32, "One saved-window job requires1..32 explicit native workers");
	}

	private void process(){
		loadLibrary();
		loadMisses();
		loadPopulation();
		loadGrade();
		loadWindows();
		loadAssociations();
		loadSequences();
		final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>();
		for(int i=0; i<Math.min(workers, selected.size()); i++){threads.add(new ProcessThread(this));}
		ThreadWaiter.startAndWait(threads);
		if(failure.get()!=null){throw new IllegalStateException("Saved-window measurement failed; no result report is accepted", failure.get());}
		for(final ProcessThread thread:threads){require(thread.success, "Each kernel worker must finish before raw attempts are joined");}
		for(final Window window:selected){
			require(window.primary!=null && window.bases==null, "Each selected window must finish one all-model primary walk and release DNA");
			primaryWork+=models.length;
			if(window.edge!=null){edgeWork+=models.length;}
		}
		reduce();
		report();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Owner Prepasses       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadLibrary(){
		models=TrnaConsensusBuilder.loadLibrary(consensus);
		names=TrnaConsensusBuilder.lastLoadedNames;
		require(models!=null && names!=null && models.length==expectedModels && names.length==models.length,
			"The diagnostic must retain every actual construction model");
		for(int i=0; i<names.length; i++){
			names[i]=RrnaResourceIO.first(names[i]);
			require(names[i].startsWith("cobalamin_RF00174_") && !modelIds.contains(names[i]), "Actual RF00174 model names must be unique");
			modelIds.put(names[i], i);
		}
	}

	private void loadMisses(){
		try(Input in=new Input(misses)){
			final byte[] header=in.reader.nextLine();
			require(header!=null && new String(header, StandardCharsets.US_ASCII).startsWith("library\toriginal_gene\tnative_truncated\tincomplete_coverage\toutcome\t"),
				"Use the reviewed full original-gene miss report");
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(!p.termEquals(label, 0) || !p.termEquals("CAPTURED_NO_GENE_ACCEPTANCE", 4)){continue;}
				final int id=p.parseInt(1);
				require(!genes.contains(id), "A library's captured/no-acceptance genes cannot repeat");
				final Gene gene=new Gene(id);
				gene.incomplete=bool(3);
				genes.put(id, gene);
				all.add(gene);
			}
		}
		require(all.size()==expectedClass2, "The saved diagnostic cannot shrink the accepted class2 population before shard association");
	}

	private void loadPopulation(){
		int count=0, incomplete=0, truncated=0, both=0;
		try(Input in=new Input(population+".original_calls.tsv.gz")){
			RrnaResourceIO.header(in.reader, PrepareNcrnaRosterBaseline.ORIGINAL_HEADER+"\tbaseline_input_id\tbaseline_input_state", population);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==31, "Original metadata must retain source coordinates and native model coverage");
				final boolean partial=bool(21), shortBoth=p.parseInt(18)-1==44 && p.parseInt(20)-p.parseInt(19)==20;
				count++;
				if(partial){incomplete++;}
				if(shortBoth){both++;}
				if(p.parseInt(17)!=0){truncated++;}
				final Gene gene=genes.get(p.parseInt(0));
				if(gene==null){continue;}
				require(!gene.population && gene.incomplete==partial, "Every rejected gene binds one unchanged original coverage record");
				gene.population=true;
				gene.shortBoth=shortBoth;
				gene.start=p.parseInt(4)-1;
				gene.stop=p.parseInt(5)-1;
				gene.reverse=p.termEquals("-", 6);
				gene.frame=p.parseString(29);
				p.setBounds(3);
				BuildNcrnaReference.field(bytes.clear(), line, p.a(), p.b());
				gene.contig=bytes.toBytes();
			}
		}
		require(count==expectedGenes && incomplete==expectedIncomplete && truncated==expectedTruncated && both==expectedShortBoth,
			"The diagnostic uses the complete accepted original, incomplete, truncated and44/20joint populations");
		for(final Gene gene:all){require(gene.population && gene.start>=0 && gene.stop>=gene.start, "No rejected original may lose its native source coordinates");}
	}

	private void loadGrade(){
		final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
		try(Input in=new Input(grade+".genes.tsv.gz")){
			RrnaResourceIO.header(in.reader, SummarizeNcrnaRosterBaseline.GENES_HEADER, grade);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				final Gene gene=genes.get(p.parseInt(0));
				if(gene==null){continue;}
				require(!gene.graded && p.termEquals(gene.frame, 1) && p.parseInt(3)-1==gene.start && p.parseInt(4)-1==gene.stop,
					"Construction-model bindings must belong to the same original gene and coordinates");
				gene.graded=true;
				if(!p.termEquals("NA", 8)){
					p.setBounds(8);
					gene.own=modelIds.get(decoder.decode(line, p.a(), p.b()));
					require(gene.own>=0, "A recorded construction model must exist in the actual complete library");
				}
			}
		}
		for(final Gene gene:all){require(gene.graded, "Every class2 gene retains its explicit construction-model binding or absence");}
	}

	private void loadWindows(){
		final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
		try(Input in=new Input(capture+".windows.tsv")){
			RrnaResourceIO.header(in.reader, CaptureNcrnaWindows.WINDOW_HEADER, capture);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==15, "Use the original captured window coordinates, including physical contig length");
				final long key=key(0, 0);
				p.setBounds(5);
				final String token=RrnaResourceIO.first(decoder.decode(line, p.a(), p.b()));
				final byte[] contig=BuildNcrnaReference.field(bytes.clear(), token).toBytes();
				final Window window=new Window(key, p.parseString(2), contig, p.parseInt(6), p.parseInt(7), p.parseInt(9), p.parseInt(10));
				require(!windows.contains(key) && window.start>=0 && window.stop>=window.start && window.stop<window.contigLength
					&& window.stop-window.start+1==p.parseInt(14) && (window.strand==0 || window.strand==1),
					"A captured event must have one valid original-oriented interval, length and strand");
				windows.put(key, window);
			}
		}
	}

	private void loadAssociations(){
		try(Input in=new Input(grade+".gene_windows.tsv.gz")){
			RrnaResourceIO.header(in.reader, SummarizeNcrnaRosterBaseline.ASSOCIATION_HEADER, grade);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				final Gene gene=genes.get(p.parseInt(0));
				if(gene==null || p.termEquals("NA", 3) || !bool(4) || !prefix(3)){continue;}
				final Window window=windows.get(key(3, capturePrefix.length));
				require(window!=null && window.frame.equals(gene.frame) && (window.strand==1)==gene.reverse && bool(8)
					&& Arrays.equals(window.contig, gene.contig), "Saved associations must keep exact event/frame/contig/strand and resolved-reference binding");
				if(gene.links.isEmpty()){local.add(gene);}
				require(!gene.links.contains(window), "An original gene/window association cannot repeat");
				gene.links.add(window);
				if(!window.selected){window.selected=true;selected.add(window);}
			}
		}
	}

	private void loadSequences(){
		Shared.TRIM_READ_DESCRIPTION=false;
		NcrnaDdlSavedData.sequences(capture+".windows.fa.gz", read->{
			bytes.clear().append(read.id);
			final Window window=windows.get(eventKey(bytes.array, 0, bytes.length()));
			require(window!=null && !window.sequenceSeen && read.bases!=null && read.length()==window.stop-window.start+1,
				"Each original captured FASTA record must bind one metadata event and retain its exact length");
			window.sequenceSeen=true;
			if(window.selected){window.bases=read.bases;}
		});
		for(final Window window:selected){require(window.bases!=null, "Every rejected gene's selected saved window needs its actual captured DNA");}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Native Kernel Workers   ----------------*/
	/*--------------------------------------------------------------*/

	private static final class ProcessThread extends Thread{
		ProcessThread(final MeasureNcrnaRejectedWindows owner_){owner=owner_;}
		@Override
		public void run(){
			try{
				for(int first=owner.next.getAndAdd(4); first<owner.selected.size(); first=owner.next.getAndAdd(4)){
					if(owner.failure.get()!=null){return;}
					for(int i=first; i<Math.min(first+4, owner.selected.size()); i++){measure(owner.selected.get(i));}
				}
				success=true;
			}catch(RuntimeException | Error invalid){owner.failure.compareAndSet(null, invalid);}
		}
		private void measure(final Window window){
			final byte[] sequence=window.bases;
			final boolean left=window.start==0, right=window.stop==window.contigLength-1;
			window.primary=new Attempts(owner.models.length);
			if(left || right){window.edge=new Attempts(owner.models.length);}
			for(int model=0; model<owner.models.length; model++){
				final float id=NcrnaAlignmentWork.selectedQuantum(owner.models[model], sequence, positions);
				window.primary.set(model, id, positions[0], positions[1], false, sequence.length);
				if(window.edge!=null){
					final EndClippedAligner.Result clipped=NcrnaAlignmentWork.contigEdges(edge, owner.models[model], sequence,
						left, right, NcrnaScavenger.EDGE_MATCH_REWARD, NcrnaScavenger.EDGE_NONMATCH_PENALTY);
					if(clipped!=null){window.edge.set(model, clipped.identity, clipped.rStart, clipped.rStop,
						clipped.qStart>0 || clipped.qStop<clipped.qLen-1, sequence.length);}
				}
			}
			window.bases=null;
		}
		final MeasureNcrnaRejectedWindows owner;
		final int[] positions=new int[4];
		final EndClippedAligner edge=new EndClippedAligner();
		boolean success;
	}

	private void reduce(){
		for(final Gene gene:local){
			for(final Window window:gene.links){
				for(int stage=0; stage<2; stage++){
					final Attempts attempts=stage==0 ? window.primary : window.edge;
					if(attempts==null){continue;}
					for(int model=0; model<models.length; model++){
						if(attempts.from[model]<0){continue;}
						final int a=window.start+attempts.from[model], b=window.start+attempts.to[model];
						final int start=window.strand==0 ? a : window.contigLength-1-b;
						final int stop=window.strand==0 ? b : window.contigLength-1-a;
						final float id=attempts.identity[model];
						final boolean overlap=overlap(gene.start, gene.stop, start, stop)>0;
						if(model==gene.own){gene.add(stage, 1, window, model, id, start, stop, attempts.genuine[model]);}
						if(!overlap){continue;}
						gene.add(stage, 0, window, model, id, start, stop, attempts.genuine[model]);
						if(model==gene.own){gene.add(stage, 2, window, model, id, start, stop, attempts.genuine[model]);}
						final boolean geometry=qualifies(gene.start, gene.stop, start, stop);
						if(!geometry || stage==1 && (!attempts.genuine[model] || id<NcrnaScavenger.EDGE_IDENTITY_FLOOR)){continue;}
						for(int cut=0; cut<3; cut++){
							if(id<CUTS[cut]){continue;}
							gene.add(stage, 3+2*cut+1, window, model, id, start, stop, attempts.genuine[model]);
							if(stop-start+1>=109){
								gene.add(stage, 3+2*cut, window, model, id, start, stop, attempts.genuine[model]);
								if(cut==2 && model==gene.own){gene.add(stage, 9, window, model, id, start, stop, attempts.genuine[model]);}
							}
						}
					}
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Reports          ----------------*/
	/*--------------------------------------------------------------*/

	private void report(){
		final long[][][] curves=new long[3][3][6];
		final LongList[][][] identities=new LongList[3][3][3];
		final LongList[][] spans=new LongList[3][3];
		final long[][][] absent=new long[3][3][3];
		for(int c=0; c<3; c++){
			for(int s=0; s<3; s++){
				for(int which=0; which<3; which++){identities[c][s][which]=new LongList();}
				spans[c][s]=new LongList();
			}
		}
		final long[] denominators=new long[3];
		try(Output out=new Output(prefix+".best.tsv.gz")){
			out.writer.println("library\tcapture\toriginal_gene\tincomplete_coverage\tmissing44left20right\tstage\tcondition\tidentity\tmodel\tspan_bases\tforward_start0\tforward_stop0\tevent\tidentity_below064\tspan_below109\tfails_reciprocal_geometry\tgenuine_clip\tconstruction_model\tedge_identity_floor_pass\twinning_attempt_stage");
			for(final Gene gene:local){
				for(int cohort=0; cohort<3; cohort++){
					if(!gene.cohort(cohort)){continue;}
					denominators[cohort]++;
					for(int stage=0; stage<3; stage++){
						for(int which=0; which<3; which++){
							final Best raw=gene.best(stage, which);
							if(raw==null){absent[cohort][stage][which]++;}
							else{
								int bin=0;
								while(bin<50 && raw.identity>=(float)((bin+1)/50.0)){bin++;}
								identities[cohort][stage][which].increment(bin);
								if(which==0){spans[cohort][stage].increment(raw.stop-raw.start+1);}
							}
						}
						for(int condition=0; condition<6; condition++){if(gene.best(stage, 3+condition)!=null){curves[cohort][stage][condition]++;}}
					}
				}
				for(int stage=0; stage<3; stage++){
					for(int condition=0; condition<CONDITIONS.length; condition++){
						final Best best=gene.best(stage, condition);
						bytes.clear().append(label).tab().append(captureId).tab().append(gene.id).tab().append(gene.incomplete).tab()
							.append(gene.shortBoth).tab().append(STAGES[stage]).tab().append(CONDITIONS[condition]).tab();
						if(best==null){bytes.append("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");}
						else{
							bytes.append(best.identity, 9).tab();
							BuildNcrnaReference.field(bytes, names[best.model]).tab();
							bytes.append(best.stop-best.start+1).tab().append(best.start).tab().append(best.stop).tab()
								.append(best.window.key).tab().append(best.identity<.64f).tab().append(best.stop-best.start+1<109).tab()
								.append(!qualifies(gene.start, gene.stop, best.start, best.stop)).tab().append(best.genuine);
						}
						bytes.tab();
						if(gene.own<0){bytes.append("NA");}else{BuildNcrnaReference.field(bytes, names[gene.own]);}
						bytes.tab();
						if(best==null || best.stage!=1){bytes.append("NA");}else{bytes.append(best.identity>=NcrnaScavenger.EDGE_IDENTITY_FLOOR);}
						bytes.tab().append(best==null ? "NA" : STAGES[best.stage]);
						out.writer.print(bytes.nl());
					}
				}
			}
		}
		try(Output curve=new Output(prefix+".curves.tsv"); Output hist=new Output(prefix+".histograms.tsv")){
			curve.writer.println("library\tcapture\tcohort\tstage\tclass2_genes\tidentity_cut\tspan_floor\tmeasured_eligible_genes_NOT_FINAL_RECALL");
			hist.writer.println("library\tcapture\tcohort\tstage\tquantity\tlower_inclusive\tupper_exclusive\tgenes");
			for(int cohort=0; cohort<3; cohort++){
				for(int stage=0; stage<3; stage++){
					for(int condition=0; condition<6; condition++){
						curve.writer.print(bytes.clear().append(label).tab().append(captureId).tab().append(COHORTS[cohort]).tab().append(STAGES[stage]).tab()
							.append(denominators[cohort]).tab().append(CUTS[condition/2], 2).tab().append(condition%2==0 ? 109 : 0).tab().append(curves[cohort][stage][condition]).nl());
					}
					for(int kind=0; kind<4; kind++){
						final int which=kind<2 ? 0 : kind-1;
						final LongList bins=kind==1 ? spans[cohort][stage] : identities[cohort][stage][which];
						long total=absent[cohort][stage][which];
						for(int bin=0; bin<bins.size; bin++){
							final long count=bins.get(bin);
							if(count==0){continue;}
							total+=count;
							bytes.clear().append(label).tab().append(captureId).tab().append(COHORTS[cohort]).tab().append(STAGES[stage]).tab().append(QUANTITIES[kind]).tab();
							if(kind!=1){bytes.append(bin/50.0, 2).tab().append(bin==50 ? "INCLUSIVE_1" : Double.toString((bin+1)/50.0));}
							else{bytes.append(bin).tab().append(bin+1);}
							hist.writer.print(bytes.tab().append(count).nl());
						}
						hist.writer.print(bytes.clear().append(label).tab().append(captureId).tab().append(COHORTS[cohort]).tab().append(STAGES[stage]).tab().append(QUANTITIES[kind])
							.append(kind==2 ? "\tMISSING_MODEL_OR_NO_ATTEMPT\tNA\t" : "\tNO_OVERLAPPING_ATTEMPT\tNA\t").append(absent[cohort][stage][which]).nl());
						require(total==denominators[cohort], "Identity and exact-base span histograms must independently conserve each shard cohort, including absent overlap");
					}
				}
			}
		}
		try(Output out=new Output(prefix+".summary.tsv")){
			out.writer.print(bytes.clear().append("metric\tcount\nfull_library_class2_genes\t").append(all.size()).append("\nshard_class2_genes\t").append(local.size())
				.append("\nselected_saved_windows\t").append(selected.size()).append("\nprimary_model_alignments\t").append(primaryWork)
				.append("\nphysical_edge_model_alignments\t").append(edgeWork).append("\nworker_count\t").append(Math.min(workers, selected.size())).nl());
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Small Helpers       ----------------*/
	/*--------------------------------------------------------------*/

	private boolean bool(final int column){
		require(p.termEquals("true", column) || p.termEquals("false", column), "Saved metadata booleans cannot be inferred");
		return p.termEquals("true", column);
	}
	private boolean prefix(final int column){
		p.setBounds(column);
		if(p.b()-p.a()<=capturePrefix.length){return false;}
		for(int i=0; i<capturePrefix.length; i++){if(p.line()[p.a()+i]!=capturePrefix[i]){return false;}}
		return true;
	}
	private long key(final int column, final int skip){p.setBounds(column);return eventKey(p.line(), p.a()+skip, p.b());}
	public static long eventKey(final byte[] id, final int from, final int to){
		require(from>=0 && from+3<to && to<=id.length && id[from]=='i', "Native event identity must be iINPUTeEVENT");
		int split=from+1;
		while(split<to && id[split]!='e'){split++;}
		require(split>from+1 && split<to-1, "Native event identity needs both original nonnegative ordinals");
		long input=0, event=0;
		for(int i=from+1; i<to; i++){
			if(i==split){continue;}
			require(id[i]>='0' && id[i]<='9', "Only decimal native event ordinals are accepted");
			if(i<split){input=input*10+id[i]-'0';}else{event=event*10+id[i]-'0';}
			require(input<=Integer.MAX_VALUE && event<=Integer.MAX_VALUE, "Native event ordinals must fit their original integer domain");
		}
		return (input<<32)|event;
	}
	private static int overlap(int a, int b, int c, int d){return Math.max(0, Math.min(b, d)-Math.max(a, c)+1);}
	private static boolean qualifies(int a, int b, int c, int d){final long n=overlap(a,b,c,d);return 2*n>=b-(long)a+1 && 2*n>=d-(long)c+1;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private static final class Attempts{
		Attempts(int count){
			identity=new float[count];
			from=new int[count];
			to=new int[count];
			genuine=new boolean[count];
			Arrays.fill(from, -1);
		}
		void set(int model, float id, int a, int b, boolean clip, int length){
			require(Float.isFinite(id) && id>=0 && id<=1 && a>=0 && b>=a && b<length, "Raw actual kernel endpoints/identity must fit their captured genomic operand");
			identity[model]=id;
			from[model]=a;
			to[model]=b;
			genuine[model]=clip;
		}
		final float[] identity;
		final int[] from, to;
		final boolean[] genuine;
	}
	private static final class Window{
		Window(long key_, String frame_, byte[] contig_, int strand_, int length_, int start_, int stop_){
			key=key_;
			frame=frame_;
			contig=contig_;
			strand=strand_;
			contigLength=length_;
			start=start_;
			stop=stop_;
		}
		final long key;
		final String frame;
		final byte[] contig;
		final int strand, contigLength, start, stop;
		boolean selected, sequenceSeen;
		byte[] bases;
		Attempts primary, edge;
	}
	private static final class Best{
		boolean better(float id, int m, Window w){return window==null || id>identity || id==identity && (m<model || m==model && w.key<window.key);}
		float identity;
		int model, start, stop, stage;
		Window window;
		boolean genuine;
	}
	private static final class Gene{
		Gene(int id_){id=id_;}
		void add(int stage, int condition, Window window, int model, float id, int start, int stop, boolean genuine){
			final int slot=stage*CONDITIONS.length+condition;
			Best best=values[slot];
			if(best==null){best=new Best();values[slot]=best;}
			if(best.better(id, model, window)){
				best.identity=id;
				best.model=model;
				best.window=window;
				best.start=start;
				best.stop=stop;
				best.genuine=genuine;
				best.stage=stage;
			}
		}
		Best best(int stage, int condition){
			if(stage<2){return values[stage*CONDITIONS.length+condition];}
			final Best a=values[condition], b=values[CONDITIONS.length+condition];
			return a==null ? b : b==null || !a.better(b.identity, b.model, b.window) ? a : b;
		}
		boolean cohort(int c){return c==0 || c==1 && incomplete || c==2 && shortBoth;}
		final int id;
		int start, stop, own=-1;
		String frame;
		byte[] contig;
		boolean incomplete, shortBoth, reverse, population, graded;
		final ArrayList<Window> links=new ArrayList<Window>();
		final Best[] values=new Best[2*CONDITIONS.length];
	}
	private String label, population, grade, capture, captureId, consensus, misses, prefix;
	private int expectedModels, expectedClass2, workers;
	private int expectedGenes=-1, expectedIncomplete=-1, expectedTruncated=-1, expectedShortBoth=-1;
	private byte[] capturePrefix;
	private byte[][] models;
	private String[] names;
	private long primaryWork, edgeWork;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder bytes=new ByteBuilder();
	private final ObjectIntMap<String> modelIds=new ObjectIntMap<String>(String.class);
	private final IntObjectMap<Gene> genes=new IntObjectMap<Gene>(16,Gene.class);
	private final LongObjectMap<Window> windows=new LongObjectMap<Window>(Window.class);
	private final ArrayList<Gene> all=new ArrayList<Gene>(), local=new ArrayList<Gene>();
	private final ArrayList<Window> selected=new ArrayList<Window>();
	private final AtomicInteger next=new AtomicInteger();
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	private static final float[] CUTS={.60f,.62f,.64f};
	private static final String[] COHORTS={"CLASS2_ALL", "CLASS2_INCOMPLETE", "CLASS2_MISSING44LEFT20RIGHT"};
	private static final String[] STAGES={"PRIMARY_RAW", "PHYSICAL_EDGE_RAW", "EITHER_RAW_NOT_FINAL_RECALL"};
	private static final String[] CONDITIONS={"BEST_POSITIVE_OVERLAP", "OWN_MODEL_ANY_OVERLAP_STATE", "OWN_MODEL_POSITIVE_OVERLAP",
		"ID060_SPAN109", "ID060_NO_SPAN", "ID062_SPAN109", "ID062_NO_SPAN", "ID064_SPAN109", "ID064_NO_SPAN", "OWN_ID064_SPAN109"};
	private static final String[] QUANTITIES={"best_overlapping_identity", "best_overlapping_span_bases", "own_model_identity_anywhere", "own_model_identity_overlapping"};
}
