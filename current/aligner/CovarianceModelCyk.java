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
import structures.IntList;
import structures.FloatList;
import structures.ListNum;
import static aligner.CovarianceModel.*;

/** Exact CYK: complete-input alignment or best target interval; global model
 * by default, optional configured local begins/EL ends. No bands, filters,
 * truncated-alignment algorithms or null3 correction.
 * Recurrence/tie order: Infernal1.1.5 cm_dpalign.c cm_CYKInsideAlign.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelCyk {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CovarianceModelCyk.class, false);final Parser p=new Parser();
		String cm=null, trace=null, coordinates=null;long maxCells=50000000;boolean local=false;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("model")){cm=value;}else if(key.equals("trace")){trace=value;}else if(key.equals("coordinates")){coordinates=value;}
			else if(key.equals("maxcells")){maxCells=Parse.parseKMG(value);}
			else if(key.equals("local")){local=Parse.parseBoolean(value);}
			else{require(p.parse(arg, key, value), "Unknown argument: "+key);}}
		require(cm!=null && p.in1!=null && p.out1!=null && trace!=null && coordinates!=null, "Bind model=, in=, out=, trace= and coordinates=");
		final String[] paths={cm,p.in1,p.out1,trace,coordinates};for(int i=0; i<paths.length; i++){for(int j=0; j<i; j++){require(!paths[i].equals(paths[j]), "Input/output paths must differ");}}
		CovarianceModelAlphabet.configureStandaloneInput();
		Shared.setThreads(1);final CovarianceModelCyk scorer=new CovarianceModelCyk(CovarianceModelParser.read(cm), maxCells, local);
		final ByteStreamWriter scores=writer(p.out1, p.overwrite), traces=writer(trace, p.overwrite), coords=writer(coordinates, p.overwrite);
		final Streamer reader=StreamerFactory.makeStreamer(FileFormat.testInput(p.in1, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 1);
		final ByteBuilder b=new ByteBuilder();final ObjectSet<String> ids=new ObjectSet<String>(String.class);int count=0, unsupported=0, impossible=0;boolean success=false;
		try{scores.println("sequence\tmodel\tlength\tstatus\tbits\tmatchedResidues\tinsertedResidues\tdeletedModelPositions\tmodelFrom1\tmodelTo1\tcells"
			+(local ? "\tlocalBeginState\telResidues\tmode" : ""));
			traces.println("sequence\tstep\tstate\ttype\ttargetFrom1\ttargetTo1\tchoice\tleftChildStep\trightChildStep");
			coords.println("sequence\ttargetPosition1\tstate\ttype\tmodelPosition1\tinsertAfterModelPosition");reader.start();
			for(ListNum<Read> batch; (batch=reader.nextList())!=null;){
				for(Read r:batch.list){final String id=first(r.id);require(ids.add(id), "Duplicate FASTA ID would make trace/score joins ambiguous: "+id);count++;final Result result;
					try{result=scorer.align(r.bases);}catch(UnsupportedSequenceException e){unsupported++;
						b.clear().append(id).tab().append(scorer.model.name).tab().append(r.length()).append("\tUNSUPPORTED_SEQUENCE\tNA\tNA\tNA\tNA\tNA\tNA\t0");
						if(local){b.append("\tNA\tNA\tlocal");}scores.print(b.nl());
						pp.outstream.println(id+": "+e.getMessage());continue;}
					if(!Float.isFinite(result.score)){impossible++;}
					b.clear().append(id).tab().append(scorer.model.name).tab().append(r.length()).tab().append(Float.isFinite(result.score) ? "OK" : "IMPOSSIBLE").tab().appendSlow(result.score)
						.tab().append(result.matched).tab().append(result.inserted).tab().append(result.deleted).tab().append(result.modelFrom).tab().append(result.modelTo).tab().append(result.cells);
					if(local){b.tab().append(result.localBeginState).tab().append(result.elResidues).append("\tlocal");}scores.print(b.nl());
					for(int n=0; n<result.states.size; n++){final int v=result.states.get(n);traces.print(b.clear().append(id).tab().append(n).tab().append(v).tab().append(scorer.model.traceStateName(v))
						.tab().append(result.from.get(n)).tab().append(result.to.get(n)).tab().append(result.choice.get(n)).tab().append(result.left.get(n)).tab().append(result.right.get(n)).nl());}
					if(Float.isFinite(result.score)){for(int i=0; i<r.length(); i++){final int v=result.stateAtBase[i];coords.print(b.clear().append(id).tab().append(i+1).tab().append(v)
							.tab().append(scorer.model.traceStateName(v)).tab().append(result.modelPosition[i]).tab().append(result.insertAnchor[i]).nl());}}
				}reader.returnList(batch);
			}require(!reader.errorState(), "FASTA streaming failed; output is incomplete");success=true;
		}finally{if(!success){reader.close();}final boolean a=scores.poisonAndWait(), z=traces.poisonAndWait(), c=coords.poisonAndWait();require(!(a || z || c), "CYK output failed");}
		pp.outstream.println("CM_CYK_PASS sequences="+count+" unsupported="+unsupported+" impossible="+impossible);Shared.closeStream(pp.outstream);
	}

	public CovarianceModelCyk(CovarianceModel model_, long maxCells_){
		this(model_, maxCells_, false);
	}
	public CovarianceModelCyk(CovarianceModel model_, long maxCells_, boolean local){
		assert(model_!=null):"The scorer requires a parsed, globally configured CM";
		require(maxCells_>0, "Positive maxcells bounds exact-matrix allocation");model=model_;maxCells=maxCells_;
		configuration=new CovarianceModelConfiguration(model, local);
		emissions=CovarianceModelAlphabet.expand(model);
		for(byte type:model.type){require(type>=D && type<=B, "Stored states must be ordinary CM states; local EL is virtual");}
	}

	/** C[0,L,L], with the entire input consumed. T and U share the RNA symbol3. */
	public Result align(byte[] sequence){return score(sequence, false, 0);}

	/** Exact best nonempty root interval, length<=maxWindow. Equal maxima retain
	 * earliest oriented end, then shortest length; every tied interval is returned.
	 * Internal DP is unchanged and still fills all intervals, without bands. */
	public Result search(byte[] sequence, int maxWindow){
		require(maxWindow>0, "Search needs an explicit positive maximum hit length to match its oracle");
		return score(sequence, true, maxWindow);
	}
	/** Best root interval with overlap >= half the shorter of candidate/hit.
	 * Candidate coordinates are one-based on the oriented window coordinate axis;
	 * they may extend beyond the window after caller boundary refinement.
	 * Both constrained and unrestricted maxima are retained, so a neighboring
	 * hit cannot silently validate a different caller candidate. */
	public Result searchOverlapping(byte[] sequence, int maxWindow, int candidateFrom, int candidateTo){
		require(maxWindow>0 && candidateFrom<=candidateTo, "Associated search requires a nonempty caller candidate and positive hit-length bound");
		return score(sequence, true, maxWindow, true, candidateFrom, candidateTo);
	}
	/** Greedy nonoverlapping local hits, each with its own score. Hard masking
	 * excludes every interval touching an extracted hit, rather than substituting
	 * ambiguous bases that a CM could still bridge. An unmasked interval's exact
	 * CYK cell depends only on its own bases, so its already-filled score is reused. */
	public Result searchGreedyOverlapping(byte[] sequence, int maxWindow, int candidateFrom, int candidateTo, float cutoff){
		require(maxWindow>0 && candidateFrom<=candidateTo && Float.isFinite(cutoff), "Greedy search needs nonempty candidate bounds and a finite model cutoff");
		return score(sequence,true,maxWindow,true,candidateFrom,candidateTo,true,cutoff);
	}

	private Result score(byte[] sequence, boolean search, int maxWindow){
		return score(sequence, search, maxWindow, false, 0, 0);
	}
	private Result score(byte[] sequence, boolean search, int maxWindow, boolean associate, int candidateFrom, int candidateTo){
		return score(sequence,search,maxWindow,associate,candidateFrom,candidateTo,false,0);
	}
	private Result score(byte[] sequence, boolean search, int maxWindow, boolean associate, int candidateFrom, int candidateTo, boolean greedy, float cutoff){
		assert(sequence!=null):"A missing sequence cannot be scored as an empty input";
		final int length=sequence.length;final long sideLong=(long)length+1, deckLong=sideLong*sideLong;
		require(deckLong<=Integer.MAX_VALUE && deckLong<=maxCells/model.states(), "Exact CYK needs states*(L+1)^2 cells; states="+model.states()+", L="+length+", maxcells="+maxCells);
		final byte[] bases=CovarianceModelAlphabet.encode(sequence);
		final int side=(int)sideLong, deck=(int)deckLong;final float[][] scores=new float[model.states()][deck];final int[][] choices=new int[model.states()][deck];
		final int[] beginState=configuration.local ? new int[deck] : null;if(beginState!=null){Arrays.fill(beginState, -1);}
		for(int v=0; v<model.states(); v++){Arrays.fill(scores[v], Float.NEGATIVE_INFINITY);Arrays.fill(choices[v], -1);}
		for(int v=model.states()-1; v>=0; v--){final byte type=model.type[v];final int delta=delta(type), rightDelta=rightDelta(type);final float[] s=scores[v];
			if(configuration.local && v==0){
				// Local root transitions are all impossible. Eligible descending
				// states accumulated their begin alternatives directly in deck0.
				for(int cell=0; cell<deck; cell++){if(beginState[cell]>=0){choices[0][cell]=LOCAL_BEGIN;}}continue;
			}
			for(int j=0; j<=length; j++){for(int d=0; d<=j; d++){final int cell=j*side+d;
				if(type==E){if(d==0){s[cell]=0;}continue;}if(d<delta){continue;}
				float best=configuration.local && Float.isFinite(configuration.endScore[v])
					? configuration.elSelf*(d-delta)+configuration.endScore[v] : Float.NEGATIVE_INFINITY;
				int choice=Float.isFinite(best) ? LOCAL_END : -1;
				if(type==B){final float[] left=scores[model.childFirst[v]], right=scores[model.childCount[v]];
					for(int k=0; k<=d; k++){final float x=left[(j-k)*side+d-k]+right[j*side+k];if(x>best){best=x;choice=k;}}}
				else{final int childCell=(j-rightDelta)*side+d-delta;
					for(int offset=0; offset<model.childCount[v]; offset++){final int child=model.childFirst[v]+offset;
						final float x=scores[child][childCell]+configuration.transitionScore[v][offset];if(x>best){best=x;choice=offset;}}
					if(delta>0){best+=emission(v, bases, j-d+1, j);}}
				s[cell]=best;choices[v][cell]=choice;
			}}
			if(configuration.local && Float.isFinite(configuration.beginScore[v])){
				for(int j=0; j<=length; j++){for(int d=0; d<=j; d++){
					final int cell=j*side+d;final float x=s[cell]+configuration.beginScore[v];
					if(x>scores[0][cell]){scores[0][cell]=x;beginState[cell]=v;}
				}}
			}
		}
		float best=scores[0][length*side+length];int from=1, to=length;
		float unrestricted=best;int unrestrictedFrom=from, unrestrictedTo=to;
		final IntList starts=new IntList(), stops=new IntList();
		if(search){best=Float.NEGATIVE_INFINITY;from=to=0;
			unrestricted=Float.NEGATIVE_INFINITY;unrestrictedFrom=unrestrictedTo=0;
			for(int j=1; j<=length; j++){for(int d=1; d<=Math.min(j, maxWindow); d++){
				final float value=scores[0][j*side+d];if(!Float.isFinite(value)){continue;}
				if(value>unrestricted){unrestricted=value;unrestrictedFrom=j-d+1;unrestrictedTo=j;}
				if(value<best || associate && !halfShorterOverlap(j-d+1, j, candidateFrom, candidateTo)){continue;}
				if(value>best){best=value;from=j-d+1;to=j;starts.clear();stops.clear();}
				starts.add(j-d+1);stops.add(j);
			}}
		}
		final GreedyHits hits=greedy ? greedyHits(scores[0],length,maxWindow,cutoff) : null;
		if(greedy){
			best=Float.NEGATIVE_INFINITY;from=to=0;starts.clear();stops.clear();
			for(int h=0;h<hits.starts.size;h++){
				if(halfShorterOverlap(hits.starts.get(h),hits.stops.get(h),candidateFrom,candidateTo)){
					if(hits.selected<0){best=hits.bits.get(h);from=hits.starts.get(h);to=hits.stops.get(h);hits.selected=h;}
					if(hits.bits.get(h)==best){starts.add(hits.starts.get(h));stops.add(hits.stops.get(h));}
				}
			}
		}
		final Result result=new Result(best, (long)deck*model.states(), length, from, to, starts, stops);
		result.greedyHits=hits;
		result.unrestrictedScore=unrestricted;result.unrestrictedFrom=unrestrictedFrom;result.unrestrictedTo=unrestrictedTo;
		if(Float.isFinite(result.score)){trace(result, scores, choices, beginState, bases, side, length);}return result;
	}
	/** Earliest end then shortest length breaks ties, matching ordinary search.
	 * The mask only removes intervals, so successive extracted scores cannot rise. */
	static GreedyHits greedyHits(float[] root,int length,int maxWindow,float cutoff){
		require(length>0 && maxWindow>0 && root.length==(long)(length+1)*(length+1) && Float.isFinite(cutoff), "Greedy selection consumes the exact root deck and finite GA threshold");
		final GreedyHits out=new GreedyHits();final boolean[] masked=new boolean[length+1];final int[] prefix=new int[length+1];final int side=length+1;
		float prior=Float.POSITIVE_INFINITY;int covered=0;
		while(covered<length){
			for(int i=1;i<=length;i++){prefix[i]=prefix[i-1]+(masked[i]?1:0);}
			float best=Float.NEGATIVE_INFINITY;int from=0,to=0;
			for(int j=1;j<=length;j++){for(int d=1;d<=Math.min(j,maxWindow);d++){
				if(prefix[j]!=prefix[j-d]){continue;}
				final float value=root[j*side+d];require(!Float.isNaN(value)&&value!=Float.POSITIVE_INFINITY,"Invalid CYK cells must not become a biological rejection");
				if(value>best){best=value;from=j-d+1;to=j;}
			}}
			if(!Float.isFinite(best)||best<cutoff){break;}
			require(from>0&&to>=from&&best<=prior,"Every greedy step must remove a real nonempty interval without increasing its score");
			out.starts.add(from);out.stops.add(to);out.bits.add(best);prior=best;
			for(int i=from;i<=to;i++){require(!masked[i],"Greedy hits cannot consume an already masked residue");masked[i]=true;covered++;}
		}
		return out;
	}
	static boolean halfShorterOverlap(int from, int to, int candidateFrom, int candidateTo){
		assert(from<=to && candidateFrom<=candidateTo):"Overlap must compare nonempty inclusive intervals";
		final long overlap=Math.max(0L, (long)Math.min(to, candidateTo)-Math.max(from, candidateFrom)+1);
		return 2*overlap>=Math.min((long)to-from+1, (long)candidateTo-candidateFrom+1);
	}

	private void trace(Result r, float[][] scores, int[][] choices, int[] beginState, byte[] bases, int side, int length){
		assert(Float.isFinite(r.score)):"Impossible parses have no valid traceback";
		final IntList stack=new IntList();push(stack, 0, r.targetTo, r.targetTo-r.targetFrom+1, -1, 0);final boolean[] matched=new boolean[model.clen+1];
		while(stack.size>0){final int sideOfParent=stack.pop(), parent=stack.pop(), d=stack.pop(), j=stack.pop(), v=stack.pop();
			final int i=j-d+1, row=r.states.size;final byte type=model.traceType(v);final int choice=type==EL ? -1 : choices[v][j*side+d];
			r.states.add(v);r.from.add(i);r.to.add(j);r.choice.add(choice);r.left.add(-1);r.right.add(-1);
			if(parent>=0){if(sideOfParent==0){r.left.set(parent, row);}else{r.right.set(parent, row);}}
			if(type==EL){
				require(configuration.local, "Virtual EL is available only in explicit local mode");
				for(int pos=i; pos<=j; pos++){
					require(r.stateAtBase[pos-1]<0, "EL cannot consume a residue already emitted by another state");
					r.stateAtBase[pos-1]=v;r.modelPosition[pos-1]=-1;r.insertAnchor[pos-1]=-1;r.elResidues++;
				}continue;
			}
			if(choice==LOCAL_BEGIN){
				require(configuration.local && v==0 && beginState!=null, "Local begin is a root-only traceback transition");
				final int entry=beginState[j*side+d];require(entry>0 && Float.isFinite(configuration.beginScore[entry]), "Selected local entry needs a configured probability");
				r.localBeginState=entry;push(stack, entry, j, d, row, 0);continue;
			}
			if(type==E){require(d==0, "An END state cannot consume residues");continue;}
			require(choice>=0 || choice==LOCAL_END && configuration.local && Float.isFinite(configuration.endScore[v]), "Finite cell is missing its traceback choice");
			if(type==B){push(stack, model.childCount[v], j, choice, row, 1);push(stack, model.childFirst[v], j-choice, d-choice, row, 0);continue;}
			final int node=model.node[v];
			if(type==MP || type==ML || type==IL){assign(r, matched, i, v, type==IL ? 0 : model.consensusLeft[node], type==IL ? model.consensusLeft[node] : -1);}
			if(type==MP || type==MR || type==IR){assign(r, matched, j, v, type==IR ? 0 : model.consensusRight[node], type==IR ? model.consensusRight[node]-1 : -1);}
			push(stack, choice==LOCAL_END ? model.states() : model.childFirst[v]+choice, j-rightDelta(type), d-delta(type), row, 0);
		}
		for(int i=1; i<=length; i++){require((r.stateAtBase[i-1]>=0)==(i>=r.targetFrom && i<=r.targetTo), "Trace must account for exactly the selected interval and leave flanks unassigned");}
		r.deleted=model.clen-r.matched;require(r.matched+r.inserted+r.elResidues==r.targetTo-r.targetFrom+1, "Trace match/insert/EL classes must conserve selected interval length");
		final float[] rescored=new float[r.states.size];
		for(int n=rescored.length-1; n>=0; n--){final int v=r.states.get(n), i=r.from.get(n), j=r.to.get(n), left=r.left.get(n), right=r.right.get(n);final byte t=model.traceType(v);
			if(t==EL){rescored[n]=configuration.elSelf*(j-i+1);continue;}
			if(t==E){rescored[n]=0;}else if(t==B){require(left>n && right>n, "B trace must retain both descendant subtrees");rescored[n]=rescored[left]+rescored[right];}
			else if(r.choice.get(n)==LOCAL_BEGIN){require(v==0 && left>n, "Local entry must join root to its selected state");rescored[n]=rescored[left]+configuration.beginScore[r.states.get(left)];}
			else{require(left>n, "Ordinary trace states must retain their chosen child");rescored[n]=rescored[left]+(r.choice.get(n)==LOCAL_END ? configuration.endScore[v] : configuration.transitionScore[v][r.choice.get(n)]);
				if(delta(t)>0){rescored[n]+=emission(v, bases, i, j);}}
			require(rescored[n]==scores[v][j*side+j-i+1], "Trace-derived score differs from its DP cell at state "+v);
		}require(rescored[0]==r.score, "Traceback must independently reconstruct the root score");
	}
	private void assign(Result r, boolean[] matched, int position, int state, int modelPosition, int anchor){
		assert(position>0 && position<=r.stateAtBase.length):"Trace emissions are one-based input positions";
		require(r.stateAtBase[position-1]<0, "Two trace states emitted the same input position");
		r.stateAtBase[position-1]=state;r.modelPosition[position-1]=modelPosition;r.insertAnchor[position-1]=anchor;
		if(modelPosition==0){require(anchor>=0 && anchor<=model.clen, "Insertion anchor leaves the consensus");r.inserted++;}
		else{require(modelPosition>0 && modelPosition<=model.clen && !matched[modelPosition], "Consensus match positions must be unique and within CLEN");matched[modelPosition]=true;r.matched++;
			if(r.modelFrom==0 || modelPosition<r.modelFrom){r.modelFrom=modelPosition;}r.modelTo=Math.max(r.modelTo, modelPosition);}
	}
	private float emission(int v, byte[] bases, int i, int j){
		assert(i>=1 && j<bases.length && i<=j):"Emitting states require a nonempty bounded input interval";
		final byte type=model.type[v];final int index=type==MP ? (bases[i]<<4)|bases[j] : (type==ML || type==IL ? bases[i] : bases[j]);return emissions[v][index];
	}
	private static void push(IntList stack, int v, int j, int d, int parent, int side){
		assert(v>=0 && d>=0 && d<=j):"Trace cells require a legal state and inclusive interval length";
		stack.add(v);stack.add(j);stack.add(d);stack.add(parent);stack.add(side);
	}
	static int delta(byte t){return t==MP ? 2 : (t==ML || t==MR || t==IL || t==IR ? 1 : 0);}
	static int rightDelta(byte t){return t==MP || t==MR || t==IR ? 1 : 0;}
	private static String first(String s){require(s!=null && !s.isEmpty(), "FASTA IDs must be nonempty");int n=0;while(n<s.length() && !Character.isWhitespace(s.charAt(n))){n++;}return s.substring(0, n);}
	private static ByteStreamWriter writer(String path, boolean overwrite){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path, FileFormat.TEXT, null, true, overwrite, false, false));w.start();return w;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	public static final class Result {
		Result(float score_, long cells_, int length, int from, int to, IntList starts, IntList stops){score=score_;cells=cells_;targetFrom=from;targetTo=to;rootStarts=starts;rootStops=stops;stateAtBase=new int[length];Arrays.fill(stateAtBase, -1);modelPosition=new int[length];insertAnchor=new int[length];Arrays.fill(insertAnchor, -1);}
		public final float score;
		public final long cells;
		/** One-based coordinates on the oriented input, inclusive;0/0 if search has no finite parse. */
		public final int targetFrom, targetTo;
		/** Unrestricted best root interval for diagnosing off-candidate winners. */
		public float unrestrictedScore;
		public int unrestrictedFrom, unrestrictedTo;
		/** Present only for explicit greedy extraction; coordinates are one-based
		 * on this input, and selected=-1 means no extracted hit associates. */
		public GreedyHits greedyHits;
		/** Exact equal-scoring best search intervals, including the selected first interval. */
		public final IntList rootStarts, rootStops;
		public final IntList states=new IntList(), from=new IntList(), to=new IntList(), choice=new IntList(), left=new IntList(), right=new IntList();
		public final int[] stateAtBase, modelPosition, insertAnchor;
		public int matched, inserted, deleted, modelFrom, modelTo;
		/** EL positions have modelPosition=-1;0 remains an ordinary insertion. */
		public int elResidues, localBeginState=-1;
	}
	public static final class GreedyHits {
		public final IntList starts=new IntList(),stops=new IntList();
		public final FloatList bits=new FloatList();
		public int selected=-1;
	}
	public static final class UnsupportedSequenceException extends IllegalArgumentException {
		UnsupportedSequenceException(String message){super(message);}
		private static final long serialVersionUID=1L;
	}
	private final CovarianceModel model;
	private final CovarianceModelConfiguration configuration;
	static final int LOCAL_BEGIN=-2, LOCAL_END=-3;
	private final float[][] emissions;
	private final long maxCells;
}
