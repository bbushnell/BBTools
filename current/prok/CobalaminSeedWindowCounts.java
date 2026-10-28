package prok;

import java.util.ArrayList;
import java.util.Locale;
import java.util.Arrays;
import dna.AminoAcid;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.FloatList;

/** Counts distinct forward seeds in the actual saved windows producing emitted
 * B12 calls. Multiple producing windows remain separate events, not extra genes.
 * Optional refs/panel/quantum inputs replay those exact pairs and compare
 * construction-style identities. Neither mode estimates counterfactual recall.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminSeedWindowCounts {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminSeedWindowCounts.class, false);
		String calls=null, stages=null, contigs=null, inputs=null, seeds=null, out=null, refs=null, panel=null;int expected=-1, quantum=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("calls")){calls=v;}else if(k.equals("stages")){stages=v;}else if(k.equals("contigs")){contigs=v;}
			else if(k.equals("inputs")){inputs=v;}else if(k.equals("seeds")){seeds=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("refs")){refs=v;}else if(k.equals("panel")){panel=v;}
			else if(k.equals("quantum")){quantum=Integer.parseInt(v);}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException(k);}
		}
		require(calls!=null && stages!=null && contigs!=null && inputs!=null && seeds!=null && out!=null && expected>0,
			"Bind one complete saved FULL_CM union call audit, its stage stream, contig roster, FASTA list and seed resource");
		Shared.setThreads(1);Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;
		require((refs==null)==(panel==null), "Paired identity diagnostics require both the exact libraries and frozen truth panel");
		require(refs==null ? quantum==-1 : quantum>=0, "Paired replay requires the caller's explicit resolved Quantum window threshold");
		try{
			final CobalaminSeedWindowCounts x=new CobalaminSeedWindowCounts();
			x.quantumThreshold=quantum;
			x.contigs(contigs);
			if(refs!=null){x.identityInputs(refs, panel);}
			x.calls(calls, expected);x.stages(stages);
			if(refs!=null){x.initialIdentities(stages);}
			x.sequences(inputs, seeds);x.write(out);
			if(refs!=null){x.writeIdentities(out);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Inputs            ----------------*/
	/*--------------------------------------------------------------*/

	private void contigs(String path){
		read(path, CobalaminPanelPrep.CONTIG_HEADER, p->{
			final String id=p.parseString(2);final int length=p.parseInt(3);
			require(length>0 && !byContig.contains(id), "Contig ownership must be unique before linking stages to FASTA");
			byContig.put(id, new Contig(p.parseString(0), id, length));
		});
	}
	private void calls(String path, int expected){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class), physical=new ObjectSet<String>(String.class);
		read(path, NcrnaReviewUnion.PANEL_CALL_HEADER, p->{
			if(!p.termEquals("FULL_CM", 0) || !p.termEquals("union", 1)){return;}
			final String id=p.parseString(2), seqid=p.parseString(4), family=p.parseString(5), model=p.parseString(9), outcome=p.parseString(10);
			final Contig c=byContig.get(seqid);final int start=p.parseInt(6), stop=p.parseInt(7), strand=p.termEquals("+", 8)?0:p.termEquals("-", 8)?1:-1;
			require(c!=null && p.termEquals(c.genome, 3) && start>0 && stop>=start && stop<=c.length && strand>=0 && ids.add(id),
				"Every emitted call must retain unique ID and exact genomic ownership/coordinates");
			require(family.equals("RF00174") || family.equals("RF01482") || family.equals("RF01689"), "Unknown B12 subtype");
			require(outcome.equals("TP") || outcome.equals("APPARENT_FP") || outcome.equals("UNSCORED"), "Unknown saved union outcome");
			require(physical.add(seqid+'\t'+start+'\t'+stop+'\t'+strand), "Exact duplicate calls require explicit allocation before window attribution");
			final Call call=new Call(id, c, family, model, outcome, start, stop, strand);c.calls.add(call);calls.add(call);
			if(compareIdentity){
				require(library.contains(model), "Selected caller model missing from paired identity library: "+model);
				if(outcome.equals("TP")){
					call.truth=truth.get(c.genome+'\t'+p.parseString(11));
					final Truth t=call.truth;
					require(t!=null && t.contig==c && t.strand==strand && t.start>0 && t.stop>=t.start,
						"A TP must retain its exact admitted truth identity, strand and core");
					final int overlap=Math.min(stop, t.stop)-Math.max(start, t.start)+1;
					require(2L*overlap>=stop-start+1L && 2L*overlap>=t.stop-t.start+1L,
						"Paired truth identity cannot transfer from a nearby unrelated locus");
				}
			}
		});
		require(calls.size()==expected, "Saved call denominator changed: "+calls.size()+" expected="+expected);
	}
	private void stages(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long line=1;
		try{
			RrnaResourceIO.header(in, NcrnaStageDiagLogSink.HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				line++;p.set(row);require(p.terms()==21, "Saved stage schema must remain21 columns");
				if(!p.termEquals("WINDOW", 0) || !p.termStartsWith("ACCEPT_", 8)){continue;}
				int family=-1;for(int i=0; i<FAMILIES.length; i++){if(p.termEquals(STAGE_FAMILIES[i], 1)){family=i;}}
				if(family<0){continue;}
				final Contig c=byContig.get(RrnaResourceIO.first(p.parseString(2)));
				require(c!=null && c.length==p.parseInt(4), "Stage contig must match the frozen panel length");
				final int strand=p.parseInt(3), pass=p.parseInt(5), a=p.parseInt(6), z=p.parseInt(7), from=p.parseInt(16), to=p.parseInt(17), raw=p.parseInt(9);
				require((strand==0 || strand==1) && (pass==1 || pass==2) && a>=0 && z>=a && z<c.length
					&& from>=0 && to>=from && to<c.length && raw>0, "Accepted stage geometry and seed occurrences must be valid");
				final int start=strand==0?from+1:c.length-to, stop=strand==0?to+1:c.length-from;
				for(Call call:c.calls){
					if(call.start==start && call.stop==stop && call.strand==strand && call.family.equals(FAMILIES[family]) && p.termEquals(call.model, 12)){
						final Window window=new Window(line, pass, a, z, raw);
						if(compareIdentity && p.termEquals("ACCEPT_IDPASS", 8)){window.direct=p.parseFloat(13);}
						call.windows.add(window);
					}
				}
			}
		}finally{require(!in.close(), "Stage input read failed");}
		for(Call c:calls){require(!c.windows.isEmpty(), "No exact accepted window/model for emitted call "+c.id+"; do not guess by proximity");}
	}
	private void sequences(String list, String seedPath){
		final LongHashSet seeds=ProkObject.loadLongKmers(seedPath, 17), matched=new LongHashSet(128);
		require(seeds!=null && seeds.size()>0, "Use the actual nonempty caller-visible forward17-mer set");
		final ByteFile in=RrnaResourceIO.open(list);final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		try{for(byte[] row; (row=in.nextLine())!=null; ){
			final String path=new String(row, java.nio.charset.StandardCharsets.UTF_8);
			require(!path.isEmpty() && paths.add(path), "The input list must preserve unique exact FASTA paths");
			RrnaResourceIO.stream(path, r->{
				final Contig c=byContig.get(RrnaResourceIO.first(r.id));
				require(c!=null && !c.seen && c.length==r.length(), "FASTA must exactly realize each frozen panel contig once");c.seen=true;
				for(Call call:c.calls){for(Window w:call.windows){
					final int raw=count(r.bases, call.strand, w.start, w.stop, seeds, matched);w.distinct=matched.size();
					require(raw==w.raw && w.distinct>0 && w.distinct<=raw, "Actual window bases must reproduce saved seed occurrences: call="+call.id+" line="+w.line+" observed="+raw+" saved="+w.raw);
					if(compareIdentity){compare(call, w, r.bases);}
				}}
			});
		}}finally{require(!in.close(), "FASTA list read failed");}
		for(Call c:calls){require(c.contig.seen, "Missing source FASTA cannot silently remove a call from the histogram");}
	}
	private void identityInputs(String refs, String panel){
		compareIdentity=true;
		for(String path:refs.split(",")){
			RrnaResourceIO.stream(path, r->{
				// loadLibrary and the GFF/stage streams retain the full FASTA header.
				final String name=r.id;
				require(!library.contains(name) && r.bases!=null && r.length()>0, "Paired replay needs unique nonempty exact models: "+name);
				library.put(name, r.bases);
			});
		}
		read(panel, CobalaminPanelPrep.HEADER, p->{
			final String key=p.parseString(1)+'\t'+p.parseString(0);final Contig c=byContig.get(p.parseString(2));
			final int a=p.parseInt(7), z=p.parseInt(8), strand=p.termEquals("+", 9)?0:p.termEquals("-", 9)?1:-1;
			require(!truth.contains(key) && c!=null && p.termEquals(c.genome, 1) && strand>=0
				&& (a==0 && z==0 || a>0 && z>=a && z<=c.length), "Truth core must retain unique genomic ownership and valid coordinates");
			truth.put(key, new Truth(c, a, z, strand));
		});
		require(library.size()>0 && truth.size()>0, "Empty libraries or truth cannot support paired identities");
	}
	/** Re-read only to bind initial MODEL observations to the selected rescue model
	 * as well as direct-pass models; WINDOW bestId can name a different model. */
	private void initialIdentities(String path){
		read(path, NcrnaStageDiagLogSink.HEADER, p->{
			if(!p.termEquals("MODEL", 0) || !p.termStartsWith("cobalamin_", 1)){return;}
			final Contig contig=byContig.get(RrnaResourceIO.first(p.parseString(2)));
			require(contig!=null && contig.length==p.parseInt(4), "MODEL observations must bind the frozen contig");
			final int strand=p.parseInt(3), pass=p.parseInt(5), a=p.parseInt(6), z=p.parseInt(7);
			for(Call c:contig.calls){if(c.strand==strand && p.termEquals(c.model, 12)){
				for(Window w:c.windows){if(w.pass==pass && w.start==a && w.stop==z){
					final float id=p.parseFloat(13);
					require(Float.isFinite(id) && id>=0 && id<=1 && (Float.isNaN(w.primary) || w.primary==id),
						"Repeated exact initial pair must have the same finite identity");w.primary=id;
				}}
			}}
		});
		for(Call c:calls){for(Window w:c.windows){
			require(Float.isFinite(w.primary) && (Float.isNaN(w.direct) || w.primary==w.direct),
				"Missing or inconsistent initial model observation for exact producing window: "+c.id);
		}}
	}
	private void compare(Call c, Window w, byte[] bases){
		final byte[] model=library.get(c.model);
		final int a=c.strand==0 ? w.start : bases.length-1-w.stop, z=c.strand==0 ? w.stop : bases.length-1-w.start;
		final byte[] window=oriented(bases, a, z, c.strand);
		final float replay=window.length>quantumThreshold ? QuantumAligner.alignStatic(model, window, positions) : ScrabbleAligner.alignStatic(window, model, null);
		require(Float.isFinite(replay) && Math.abs(replay-w.primary)<=0.000001f,
			"Saved generic primary identity must replay on the exact window/model; do not substitute another aligner or geometry: "+c.id);
		w.scrabble=ScrabbleAligner.alignStatic(window, model, null);
		if(c.truth!=null){
			final byte[] core=oriented(bases, c.truth.start-1, c.truth.stop-1, c.strand);
			w.member=ScrabbleAligner.alignStatic(core, model, null);
		}
		assert(Float.isFinite(w.scrabble) && (c.truth==null || Float.isFinite(w.member))) : "Paired alignment metrics must remain finite before reporting deltas";
	}
	public static byte[] oriented(byte[] bases, int start, int stop, int strand){
		require(bases!=null && start>=0 && stop>=start && stop<bases.length && (strand==0 || strand==1), "Extract only valid inclusive genomic intervals");
		final byte[] copy=Arrays.copyOfRange(bases, start, stop+1);shared.Tools.toUpperCase(copy);
		if(strand==1){AminoAcid.reverseComplementBasesInPlace(copy);}return copy;
	}
	private void writeIdentities(String out){
		RrnaResourceIO.fresh(out, ".identities.tsv", ".identity_summary.tsv");
		final ByteStreamWriter rows=RrnaResourceIO.writer(out+".identities.tsv"), summary=RrnaResourceIO.writer(out+".identity_summary.tsv");
		final ByteBuilder b=new ByteBuilder();final FloatList[][] deltas=new FloatList[4][2];
		for(int f=0; f<4; f++){for(int d=0; d<2; d++){deltas[f][d]=new FloatList();}}
		try{
			rows.println("call_id\tfamily\tmodel\toutcome\tstage_line\twindow_length\ttruth_core_length\tcaller_primary\tScrabble_same_window\tScrabble_truth_core\tprimary_minus_window_Scrabble\tprimary_minus_core_Scrabble");
			for(Call c:calls){int f=-1;for(int i=0; i<3; i++){if(c.family.equals(FAMILIES[i])){f=i;}}
				assert(f>=0) : "Calls were validated against the three B12 subfamilies before pairing";
				for(Window w:c.windows){
					final float windowDelta=w.primary-w.scrabble;
					b.clear().append(c.id).tab().append(c.family).tab().append(c.model).tab().append(c.outcome).tab().append(w.line)
						.tab().append(w.stop-w.start+1).tab().append(c.truth==null ? 0 : c.truth.stop-c.truth.start+1)
						.tab().append(w.primary, 9).tab().append(w.scrabble, 9).tab();
					if(c.truth==null){b.append("NA");}else{b.append(w.member, 9);}
					b.tab().append(windowDelta, 9).tab();
					if(c.truth==null){b.append("NA");}else{
						final float coreDelta=w.primary-w.member;b.append(coreDelta, 9);
						deltas[f][1].add(coreDelta);deltas[3][1].add(coreDelta);
					}
					rows.print(b.nl());deltas[f][0].add(windowDelta);deltas[3][0].add(windowDelta);
				}
			}
			summary.println("family\tcomparison\tpaired_window_events\tmin\tp5\tmedian\tp95\tmax\tmean\tscope");
			for(int f=0; f<4; f++){for(int d=0; d<2; d++){
				final FloatList values=deltas[f][d];values.sort();double sum=0;for(int i=0; i<values.size; i++){sum+=values.array[i];}
				b.clear().append(f==3 ? "union" : FAMILIES[f]).tab().append(d==0 ? "primary_minus_same_window_Scrabble" : "primary_minus_truth_core_Scrabble").tab().append(values.size);
				for(int q:new int[]{0,5,50,95,100}){b.tab();if(values.size==0){b.append("NA");}else{b.append(values.array[(int)((values.size-1L)*q/100)], 9);}}
				b.tab();if(values.size==0){b.append("NA");}else{b.append(sum/values.size, 9);}
				summary.print(b.append(d==0 ? "\tEMITTED_CALLS_ONLY_NOT_ALL_ATTEMPTED_MODELS\n" : "\tMATCHED_TP_CM_CORES_SELECTED_MODELS_NOT_ORIGINAL_TRAINING_ASSIGNMENTS\n"));
			}}
		}finally{RrnaResourceIO.finish(rows, summary);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(String out){
		RrnaResourceIO.fresh(out, ".windows.tsv", ".calls.tsv", ".histogram.tsv");
		final ByteStreamWriter windows=RrnaResourceIO.writer(out+".windows.tsv"), summary=RrnaResourceIO.writer(out+".calls.tsv"), hist=RrnaResourceIO.writer(out+".histogram.tsv");
		final ByteBuilder b=new ByteBuilder();final long[][][] counts=new long[4][3][6];int events=0;
		try{
			windows.println("call_index\tgenome\tseqid\temitted_family\toutcome\tstrand\tstage_line\tpass\twindow_start0_oriented\twindow_stop0_oriented\tseed_occurrences\tdistinct_seed_keys");
			summary.println("call_index\temitted_family\toutcome\tproducing_windows\tminimum_distinct_keys\tmaximum_distinct_keys");
			for(Call c:calls){
				int minimum=Integer.MAX_VALUE, maximum=0, family=-1;for(int i=0; i<FAMILIES.length; i++){if(c.family.equals(FAMILIES[i])){family=i;}}
				final int outcome=c.outcome.equals("TP")?0:c.outcome.equals("APPARENT_FP")?1:2;
				for(Window w:c.windows){
					require(w.distinct>0, "Only measured positive-seed windows enter the histogram");events++;
					minimum=Math.min(minimum, w.distinct);maximum=Math.max(maximum, w.distinct);final int bucket=Math.min(5, w.distinct);
					counts[family][outcome][bucket]++;counts[3][outcome][bucket]++;
					windows.print(b.clear().append(c.id).tab().append(c.contig.genome).tab().append(c.contig.id).tab().append(c.family).tab().append(c.outcome)
						.tab().append(c.strand==0?'+':'-').tab().append(w.line).tab().append(w.pass).tab().append(w.start).tab().append(w.stop).tab().append(w.raw).tab().append(w.distinct).nl());
				}
				summary.print(b.clear().append(c.id).tab().append(c.family).tab().append(c.outcome).tab().append(c.windows.size()).tab().append(minimum).tab().append(maximum).nl());
			}
			hist.println("scope\toutcome\tdistinct_seed_bucket\tproducing_window_events");long total=0;
			for(int f=0; f<4; f++){for(int o=0; o<3; o++){for(int k=1; k<=5; k++){
				final long n=counts[f][o][k];if(f==3){total+=n;}
				hist.print(b.clear().append(f==3?"union":FAMILIES[f]).tab().append(OUTCOMES[o]).tab().append(k==5?"5+":Integer.toString(k)).tab().append(n).nl());
			}}}
			require(total==events, "Histogram must conserve all producing-window events, without reporting them as additional genes");
		}finally{RrnaResourceIO.finish(windows, summary, hist);}
		System.err.println("COBALAMIN_SEED_WINDOW_COUNTS_PASS calls="+calls.size()+" producingWindows="+events+"; diagnostic only, not two-seed counterfactual accuracy");
	}
	/** Counts in the saved strand-oriented inclusive window without copying it. */
	public static int count(byte[] bases, int strand, int start, int stop, LongHashSet seeds, LongHashSet matched){
		require(bases!=null && (strand==0 || strand==1) && start>=0 && stop>=start && stop<bases.length && seeds!=null && matched!=null,
			"Seed recount requires the exact genomic sequence and valid oriented window");
		matched.clear();long word=0;int valid=0, raw=0;
		for(int i=start; i<=stop; i++){
			int code=RrnaResourceIO.base(bases[strand==0?i:bases.length-1-i]);
			if(code<0){word=0;valid=0;continue;}if(strand==1){code=3-code;}
			word=((word<<2)|code)&MASK;
			if(++valid>=17 && seeds.contains(word)){raw++;matched.add(word);}
		}
		assert(matched.size()<=raw):"Each distinct forward key needs at least one in-window occurrence";return raw;
	}
	private static void read(String path, String header, java.util.function.Consumer<LineParser1> action){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');p.set(header.getBytes(java.nio.charset.StandardCharsets.UTF_8));final int terms=p.terms();
		try{RrnaResourceIO.header(in, header, path);for(byte[] row; (row=in.nextLine())!=null; ){p.set(row);require(p.terms()==terms, "Input table schema changed: "+path);action.accept(p);}}
		finally{require(!in.close(), "Input table read failed: "+path);}
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Contig {
		Contig(String g, String i, int n){genome=g;id=i;length=n;}
		final String genome, id;final int length;boolean seen;final ArrayList<Call> calls=new ArrayList<Call>();
	}
	private static final class Call {
		Call(String i, Contig c, String f, String m, String o, int a, int z, int s){id=i;contig=c;family=f;model=m;outcome=o;start=a;stop=z;strand=s;}
		final String id, family, model, outcome;final Contig contig;final int start, stop, strand;final ArrayList<Window> windows=new ArrayList<Window>();Truth truth;
	}
	private static final class Window {
		Window(long n, int p, int a, int z, int r){line=n;pass=p;start=a;stop=z;raw=r;}
		final long line;final int pass, start, stop, raw;int distinct;
		float primary=Float.NaN, direct=Float.NaN, scrabble=Float.NaN, member=Float.NaN;
	}
	private static final class Truth {
		Truth(Contig c, int a, int z, int s){contig=c;start=a;stop=z;strand=s;}
		final Contig contig;final int start, stop, strand;
	}
	private boolean compareIdentity=false;
	private int quantumThreshold=-1;
	private final int[] positions=new int[2];
	private final ObjectMap<String,byte[]> library=new ObjectMap<String,byte[]>(String.class, byte[].class);
	private final ObjectMap<String,Truth> truth=new ObjectMap<String,Truth>(String.class, Truth.class);
	private final ObjectMap<String,Contig> byContig=new ObjectMap<String,Contig>(String.class, Contig.class);
	private final ArrayList<Call> calls=new ArrayList<Call>();
	private static final long MASK=(1L<<34)-1;
	private static final String[] FAMILIES={"RF00174", "RF01482", "RF01689"};
	private static final String[] STAGE_FAMILIES={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
	private static final String[] OUTCOMES={"TP", "APPARENT_FP", "UNSCORED"};
}
