package aligner;

import java.util.Arrays;
import java.util.Locale;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import static aligner.CovarianceModel.*;

/** Exact full-input CYK with triangular score decks and parent-count reuse.
 * Recurrence and float operation order match CovarianceModelCyk; no traceback,
 * bands, truncation, search or null3 correction. Optional local model begins
 * and EL ends use the same immutable configuration as the dense scorer.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelScoreOnly {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CovarianceModelScoreOnly.class, false);final Parser p=new Parser();
		String cm=null;long maxCells=50000000;boolean local=false;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("model")){cm=value;}else if(key.equals("maxcells")){maxCells=Parse.parseKMG(value);}
			else if(key.equals("local")){local=Parse.parseBoolean(value);}
			else{require(p.parse(arg, key, value), "Unknown argument: "+key);}}
		require(cm!=null && p.in1!=null && p.out1!=null, "Bind model=, in= and out=");
		require(!cm.equals(p.in1) && !cm.equals(p.out1) && !p.in1.equals(p.out1), "Input and output paths must differ");
		CovarianceModelAlphabet.configureStandaloneInput();
		Shared.setThreads(1);final CovarianceModelScoreOnly scorer=new CovarianceModelScoreOnly(CovarianceModelParser.read(cm), maxCells, local);
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(p.out1, FileFormat.TEXT, null, true, p.overwrite, false, false));
		final Streamer reader=StreamerFactory.makeStreamer(FileFormat.testInput(p.in1, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 1);
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);final ByteBuilder b=new ByteBuilder();
		int count=0, unsupported=0, impossible=0;boolean success=false;out.start();
		try{out.println("sequence\tmodel\tlength\tstatus\tbits\tvisitedCells\tdeckCells\tallocatedDecks\tpeakScoreCells\tseconds"+(local ? "\tmode" : ""));reader.start();
			for(ListNum<Read> batch; (batch=reader.nextList())!=null;){for(Read r:batch.list){
				final String id=first(r.id);require(ids.add(id), "Duplicate FASTA ID makes score joins ambiguous: "+id);count++;
				final Result result;final long start=System.nanoTime();
				try{result=scorer.align(r.bases);}catch(CovarianceModelCyk.UnsupportedSequenceException e){unsupported++;
					b.clear().append(id).tab().append(scorer.model.name).tab().append(r.length()).append("\tUNSUPPORTED_SEQUENCE\tNA\t0\t0\t0\t0\tNA");
					if(local){b.append("\tlocal");}out.print(b.nl());
					pp.outstream.println(id+": "+e.getMessage());continue;}
				final double seconds=(System.nanoTime()-start)*1e-9;if(!Float.isFinite(result.score)){impossible++;}
				b.clear().append(id).tab().append(scorer.model.name).tab().append(r.length()).tab().append(Float.isFinite(result.score) ? "OK" : "IMPOSSIBLE")
					.tab().appendSlow(result.score).tab().append(result.visitedCells).tab().append(result.deckCells).tab().append(result.allocatedDecks)
					.tab().append(result.peakScoreCells).tab().append(seconds, 6);
				if(local){b.append("\tlocal");}out.print(b.nl());
			}reader.returnList(batch);}
			require(!reader.errorState(), "FASTA streaming failed; score output is incomplete");success=true;
		}finally{if(!success){reader.close();}require(!out.poisonAndWait(), "Score-only output failed");}
		pp.outstream.println("CM_SCORE_ONLY_PASS sequences="+count+" unsupported="+unsupported+" impossible="+impossible);Shared.closeStream(pp.outstream);
	}

	public CovarianceModelScoreOnly(CovarianceModel model_, long maxCells_){
		this(model_, maxCells_, false);
	}
	public CovarianceModelScoreOnly(CovarianceModel model_, long maxCells_, boolean local){
		require(model_!=null && maxCells_>0, "Score-only requires a parsed model and a positive resident-cell limit");
		model=model_;maxCells=maxCells_;parents=new int[model.states()];
		configuration=new CovarianceModelConfiguration(model, local);
		emissionScores=CovarianceModelAlphabet.expand(model);
		for(int v=0; v<model.states(); v++){
			final byte type=model.type[v];require(type>=D && type<=B, "Stored CM states exclude the virtual local-end state");
			if(type==B){addParent(v, model.childFirst[v]);addParent(v, model.childCount[v]);}
			else if(type==E){require(model.childCount[v]==0, "END states cannot consume child decks");}
			else{require(model.childCount[v]>0 && model.transitionScore[v].length==model.childCount[v], "Transition scores must cover the child range");
				for(int c=0; c<model.childCount[v]; c++){addParent(v, model.childFirst[v]+c);}}
		}
		// Derive lifetimes from the edges actually read by the recurrence, not header metadata.
		final int[] pending=parents.clone();int live=0, peak=0;
		for(int v=model.states()-1; v>=0; v--){live++;peak=Math.max(peak, live);
			if(model.type[v]==B){live-=consume(v, model.childFirst[v], pending);live-=consume(v, model.childCount[v], pending);}
			else if(model.type[v]!=E){for(int c=0; c<model.childCount[v]; c++){live-=consume(v, model.childFirst[v]+c, pending);}}
			if(v!=0 && pending[v]==0){live--;}
		}
		require(live==1, "Deck lifetime plan must retain only the root at completion");peakDecks=peak;
	}

	/** Score the entire IUPAC input; maxcells bounds resident float cells, not work. */
	public Result align(byte[] sequence){
		require(sequence!=null, "Missing sequence cannot be interpreted as empty input");
		final int length=sequence.length;final long cells=triangle(length);
		require(cells<=Integer.MAX_VALUE-8 && cells<=maxCells/peakDecks,
			"Triangular score decks exceed resident-cell limit: L="+length+", decks="+peakDecks+", maxcells="+maxCells);
		final byte[] bases=CovarianceModelAlphabet.encode(sequence);final int deck=(int)cells;final int[] rows=new int[length+1];
		for(int j=0; j<=length; j++){rows[j]=(int)((long)j*(j+1)/2);}
		final float[][] scores=new float[model.states()][], pool=new float[peakDecks][];
		final int[] pending=parents.clone();int available=0, allocated=0;
		float localBest=Float.NEGATIVE_INFINITY;
		for(int v=model.states()-1; v>=0; v--){
			final float[] s;if(available>0){s=pool[--available];pool[available]=null;}else{s=new float[deck];allocated++;}
			Arrays.fill(s, Float.NEGATIVE_INFINITY);scores[v]=s;
			final byte type=model.type[v];final int delta=CovarianceModelCyk.delta(type), rightDelta=CovarianceModelCyk.rightDelta(type);
			if(type==E){for(int j=0; j<=length; j++){s[rows[j]]=0;}}
			else if(type==B){final float[] left=scores[model.childFirst[v]], right=scores[model.childCount[v]];
				assert(left!=null && right!=null):"Pending bifurcation sides must survive until their parent is scored, state="+v;
				for(int j=0; j<=length; j++){final int row=rows[j];for(int d=0; d<=j; d++){
					float best=Float.NEGATIVE_INFINITY;
					for(int k=0; k<=d; k++){final float x=left[rows[j-k]+d-k]+right[row+k];if(x>best){best=x;}}
					s[row+d]=best;
				}}
			}else{final int first=model.childFirst[v], count=model.childCount[v];final float[] transitions=configuration.transitionScore[v], emissions=emissionScores[v];
				for(int j=delta; j<=length; j++){final int row=rows[j], childRow=rows[j-rightDelta];for(int d=delta; d<=j; d++){
					float best=configuration.local && Float.isFinite(configuration.endScore[v])
						? configuration.elSelf*(d-delta)+configuration.endScore[v] : Float.NEGATIVE_INFINITY;
					final int childCell=childRow+d-delta;
					for(int c=0; c<count; c++){
						assert(scores[first+c]!=null):"A child deck was released before its last parent consumed it, state="+v+", child="+(first+c);
						final float x=scores[first+c][childCell]+transitions[c];if(x>best){best=x;}}
					if(delta>0){final int i=j-d+1, index=type==MP ? (bases[i]<<4)|bases[j] : (type==ML || type==IL ? bases[i] : bases[j]);best+=emissions[index];}
					s[row+d]=best;
				}}
			}
			// Capture the full-input local-entry alternative before this deck can
			// be recycled; no extra retained matrix or traceback is required.
			if(configuration.local && Float.isFinite(configuration.beginScore[v])){
				localBest=Math.max(localBest, s[rows[length]+length]+configuration.beginScore[v]);
			}
			if(type==B){if(consume(v, model.childFirst[v], pending)>0){pool[available++]=scores[model.childFirst[v]];scores[model.childFirst[v]]=null;}
				if(consume(v, model.childCount[v], pending)>0){pool[available++]=scores[model.childCount[v]];scores[model.childCount[v]]=null;}}
			else if(type!=E){for(int c=0; c<model.childCount[v]; c++){final int child=model.childFirst[v]+c;
				if(consume(v, child, pending)>0){pool[available++]=scores[child];scores[child]=null;}}}
			if(v!=0 && pending[v]==0){pool[available++]=s;scores[v]=null;}
		}
		require(allocated==peakDecks && available==allocated-1, "Actual deck reuse disagrees with the topology lifetime plan");
		return new Result(configuration.local ? localBest : scores[0][rows[length]+length], cells*model.states(), deck, allocated);
	}

	private void addParent(int v, int child){
		require(child>=v && child<model.states(), "Descending-state evaluation requires forward child indices: "+v+" -> "+child);
		if(child==v){require(model.type[v]==IL || model.type[v]==IR, "Only consuming insertion states can depend on their own earlier cells");}
		else{parents[child]++;}
	}
	private static int consume(int v, int child, int[] pending){
		if(child==v){return 0;}require(pending[child]>0, "Every child edge must be consumed once before its deck can be reused: "+v+" -> "+child);
		return --pending[child]==0 ? 1 : 0;
	}
	static long triangle(int length){
		require(length>=0, "Triangular intervals require a nonnegative sequence length");return ((long)length+1)*(length+2L)/2;
	}
	private static String first(String s){require(s!=null && !s.isEmpty(), "FASTA IDs must be nonempty");int n=0;while(n<s.length() && !Character.isWhitespace(s.charAt(n))){n++;}return s.substring(0, n);}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	public static final class Result {
		Result(float score_, long visited_, int deck_, int allocated_){score=score_;visitedCells=visited_;deckCells=deck_;allocatedDecks=allocated_;peakScoreCells=(long)deck_*allocated_;}
		public final float score;
		/** Logical triangular cells, including initialized impossible cells; not B split evaluations. */
		public final long visitedCells, peakScoreCells;
		public final int deckCells, allocatedDecks;
	}
	private final CovarianceModel model;
	private final CovarianceModelConfiguration configuration;
	private final float[][] emissionScores;
	private final long maxCells;
	private final int[] parents;
	final int peakDecks;
}
