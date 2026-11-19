package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import cardinality.DynamicDemiLog;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData.Edge;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Diagnostic replay of frozen ncRNA windows after the shared N-boundary repair.
 * The deliberately preserved old loop is used only to reproduce saved controls.
 * It is not an alternative production encoder. Three states separate query
 * changes from member-profile changes: old/old, fixed-query/old-profile, fixed/fixed.
 * @author Ganyu
 */
public final class CompareNcrnaHashBoundary{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, CompareNcrnaHashBoundary.class, false);
		final int exponent=DynamicDemiLog.exponentBits();
		try{
			DynamicDemiLog.setExponent(4);
			final Config c=new Config();
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("panels")){c.panels=b;}
				else if(a.equals("selection")){c.selection=b;}
				else if(a.equals("baseline")){c.baseline=b;}
				else if(a.equals("members")){c.members=b;}
				else if(a.equals("ledger")){c.ledger=b;}
				else if(a.equals("assignments")){c.assignments=b;}
				else if(a.equals("census")){c.census=b;}
				else if(a.equals("summary")){c.summary=b;}
				else if(a.equals("ref")){c.references=b;}
				else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
				else if(a.equals("expectedgroups")){c.groups=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown boundary replay option: "+arg);}
			}
			require(parser.in1==null && parser.in2==null && parser.out1!=null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use complete saved bindings and out=fresh_prefix, never genomes or sampling");
			c.output=parser.out1;
			c.threads=Shared.threads();
			run(c);
			pp.outstream.println("NCRNA_HASH_BOUNDARY_SAVED_COMPARISON_COMPLETE");
		}finally{
			DynamicDemiLog.setExponent(exponent);
			Shared.closeStream(pp.outstream);
		}
	}

	public static void run(final Config c){
		require(c!=null && c.panels!=null && c.selection!=null && c.baseline!=null && c.members!=null
			&& c.ledger!=null && c.assignments!=null && c.census!=null && c.summary!=null
			&& c.references!=null && c.output!=null && c.models>0 && c.groups>0 && c.threads>0,
			"Replay requires complete frozen cache, old hybrid result, member bindings and positive sizes");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions && DynamicDemiLog.exponentBits()==4 && !DynamicDemiLog.blacklistExists(),
			"Replay needs -ea, native exponent4 and no blacklist");
		RrnaResourceIO.fresh(c.output, ".profiles.tsv", ".changes.tsv.gz", ".loci.tsv.gz", ".summary.tsv");
		final NcrnaDdlQueryKeys fixedGuard=new NcrnaDdlQueryKeys(32, 4);
		fixedGuard.build("AAANAAA".getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		require(fixedGuard.keyCount()==0, "Run this comparison only on the fixed encoder, never mislabel an old build as fixed");
		Shared.TRIM_READ_DESCRIPTION=false;
		new CompareNcrnaHashBoundary(c).process();
	}

	private CompareNcrnaHashBoundary(final Config c_){c=c_;}

	/** Literal diagnostic used by the test: one member per model, same three raw-score
	 * paths as the real replay. The last model is the independently specified target. */
	public static int[] literalRanks(final byte[][] members, final byte[] bases, final int k, final int width){
		require(members!=null && members.length>0 && bases!=null, "Literal comparison requires members and one query");
		final NcrnaDdlProfile[][] profiles=new NcrnaDdlProfile[2][members.length];
		final DynamicDemiLog sketch=DynamicDemiLog.create(width, k, 12345, 0, false, false);
		for(int model=0; model<members.length; model++){
			for(int state=0; state<2; state++){
				sketch.clear();
				if(state==0){legacy(members[model], k, true, false, sketch, null);}
				else{sketch.hash(members[model], null);}
				final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(width);
				votes.add(sketch.maxArray());
				profiles[state][model]=votes.finish();
			}
		}
		final LegacyRanks old=new LegacyRanks(profiles[0], k, true, false);
		final NcrnaDdlRanks query=new NcrnaDdlRanks(profiles[0], k);
		final NcrnaDdlRanks fixed=new NcrnaDdlRanks(profiles[1], k);
		old.score(bases);
		query.score(bases);
		fixed.score(bases);
		final boolean[] acceptable=new boolean[members.length];
		acceptable[members.length-1]=true;
		final int[] tuple=new int[4];
		final int[] result=new int[3];
		old.best(NcrnaDdlRanks.ALL_KEYS, acceptable, tuple);
		result[0]=tuple[0];
		query.best(NcrnaDdlRanks.ALL_KEYS, acceptable, tuple);
		result[1]=tuple[0];
		fixed.best(NcrnaDdlRanks.ALL_KEYS, acceptable, tuple);
		result[2]=tuple[0];
		return result;
	}

	private void process(){
		data=NcrnaDdlSavedData.load(c.panels, c.selection, c.models, c.groups);
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(c.references, read->{
			names.add(read.id);
			sequences.add(read.bases);
		});
		require(names.size()==c.models, "References must preserve the complete cached model order");
		final String[] labels=names.toArray(new String[names.size()]);
		final byte[][] labelBytes=new byte[labels.length][];
		for(int model=0; model<labels.length; model++){
			labelBytes[model]=labels[model].getBytes(java.nio.charset.StandardCharsets.UTF_8);
		}
		final NcrnaDdlMembers members=NcrnaDdlMembers.load(c.members, c.ledger, c.assignments,
			c.census, c.summary, labels, sequences.toArray(new byte[sequences.size()][]));
		data.readNativeCounts(labels);
		baseline();
		profiles(members, labelBytes);
		windows=data.windows.toArray(new Window[data.windows.size()]);
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		try(Output changes=new Output(c.output+".changes.tsv.gz")){
			changes.writer.println("panel\tinput_slot\tinput_id\tevent\tparent\tseen_stratum\tquery_Ns\tmethod\told_rank\told_tie_first\told_tie_last\told_model\tquery_fixed_rank\tquery_fixed_tie_first\tquery_fixed_tie_last\tquery_fixed_model\tall_fixed_rank\tall_fixed_tie_first\tall_fixed_tie_last\tall_fixed_model");
			for(int i=0; i<Math.min(c.threads, Math.max(1, windows.length)); i++){workers.add(new Worker(changes));}
			ThreadWaiter.startAndWait(workers);
			if(failure.get()!=null){throw new IllegalStateException("Saved ambiguity comparison failed before publication", failure.get());}
		}
		report(workers);
	}

	/** Original accepted hybrid locus table; all three old minima remain strict controls. */
	private void baseline(){
		oldMin=new int[data.parents.size()*3];
		final boolean[] seen=new boolean[data.parents.size()];
		final ObjectMap<String, Parent> parents=new ObjectMap<String, Parent>(String.class, Parent.class);
		for(final Parent parent:data.parents){parents.put(parent.input.panel.name+'\t'+new String(parent.idBytes, java.nio.charset.StandardCharsets.UTF_8), parent);}
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(c.baseline)){
			NcrnaDdlSavedData.header(in, "panel\tinput_slot\tinput_id\tparent\tseen_stratum\tnative_rank\tsketch_rank\thybrid_rank\tsource_disposition");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==9, "Old hybrid locus table must retain its exact schema");
				final String panel=NcrnaDdlSavedData.text(p, 0);
				final String parentId=NcrnaDdlSavedData.text(p, 3);
				final Parent parent=parents.get(panel+'\t'+parentId);
				require(parent!=null && !seen[parent.index] && p.parseInt(1)==parent.input.slot
					&& NcrnaDdlSavedData.text(p, 2).equals(parent.input.id) && p.termEquals(NcrnaDdlSavedData.STRATA[parent.stratum], 4)
					&& p.termEquals(parent.reason, 8), "Old minima must bind every original selected parent exactly once; panel="
					+panel+", parent="+parentId+", input="+NcrnaDdlSavedData.text(p, 2));
				seen[parent.index]=true;
				for(int m=0; m<3; m++){
					final int rank=p.parseInt(5+m);
					require(rank==-1 || rank>0 && rank<=c.models, "Old minima must be explicit ranks or unrankable");
					oldMin[parent.index*3+m]=rank;
				}
				require(oldMin[parent.index*3]==parent.nativeRank, "Old hybrid cannot replace the original native minimum");
			}
		}
		for(final boolean present:seen){require(present, "Old hybrid must include zero-support and unrankable targets");}
	}

	private void profiles(final NcrnaDdlMembers members, final byte[][] labels){
		try(Output out=new Output(c.output+".profiles.tsv")){
			out.writer.println("geometry\tk\tbuckets\tcanonical\tmiddle_mask\tmodel_index\tmodel\tmembers\tN_containing_members\tchanged_key_buckets");
			final ByteBuilder row=new ByteBuilder();
			for(int g=0; g<9; g++){
				final int k=k(g);
				final int width=width(g);
				fresh[g]=NcrnaDdlModelProfiles.build(members, k, width, g>0, g==0);
				old[g]=new NcrnaDdlProfile[c.models];
				final DynamicDemiLog sketch=DynamicDemiLog.create(width, k, 12345, 0, false, false);
				for(int model=0; model<c.models; model++){
					final NcrnaDdlProfile.Votes votes=new NcrnaDdlProfile.Votes(width);
					long nMembers=0;
					for(final byte[] bases:members.members(model)){
						if(countN(bases)>0){nMembers++;}
						sketch.clear();
						legacy(bases, k, g>0, g==0, sketch, null);
						votes.add(sketch.maxArray());
					}
					old[g][model]=votes.finish();
					int changed=0;
					for(int b=0; b<width; b++){if(old[g][model].keys()[b]!=fresh[g][model].keys()[b]){changed++;}}
					require(nMembers>0 || changed==0, "N-free member profiles must be byte-identical");
					out.writer.print(profileRow(row, g, model, labels[model], members.members(model).length, nMembers, changed));
				}
			}
		}
	}

	/** Actual profile report formatter; retained UTF-8 model bytes never pass through
	 * ByteBuilder.append(String), which narrows Java characters to bytes. */
	public static ByteBuilder profileRow(final ByteBuilder row, final int geometry, final int model,
			final byte[] label, final int members, final long nMembers, final int changed){
		require(row!=null && geometry>=0 && geometry<9 && model>=0 && label!=null && label.length>0
			&& members>0 && nMembers>=0 && nMembers<=members && changed>=0 && changed<=width(geometry),
			"Profile report requires one valid geometry/model and conserved member/bucket counts");
		return row.clear().append(geometry).tab().append(k(geometry)).tab().append(width(geometry)).tab()
			.append(geometry>0 ? 1 : 0).tab().append(geometry==0 ? 1 : 0).tab().append(model).tab().append(label).tab()
			.append(members).tab().append(nMembers).tab().append(changed).nl();
	}

	private final class Worker extends Thread{
		Worker(final Output output){
			changes=output;
			Arrays.fill(minima, -1);
			for(int g=0; g<9; g++){
				legacyRanks[g]=new LegacyRanks(old[g], k(g), g>0, g==0);
				queryRanks[g]=new NcrnaDdlRanks(old[g], k(g), g>0, g==0);
				fixedRanks[g]=new NcrnaDdlRanks(fresh[g], k(g), g>0, g==0);
			}
			for(int state=0; state<3; state++){order[state]=new NcrnaDdlInterleave(c.models);}
		}
		@Override public void run(){
			try{
				for(int index=next.getAndIncrement(); index<windows.length && failure.get()==null; index=next.getAndIncrement()){window(windows[index]);}
				flush();
				success=true;
			}catch(Throwable problem){failure.compareAndSet(null, problem);}
		}
		private void window(final Window window){
			final int ns=countN(window.bases);
			windowCounts[ns==0 ? 0 : 1]++;
			boolean useful=false;
			for(final Edge edge:window.edges){useful|=edge.unique>0;}
			if(useful){
				for(int g=0; g<9; g++){
					legacyRanks[g].score(window.bases);
					queryRanks[g].score(window.bases);
					fixedRanks[g].score(window.bases);
				}
				for(int state=0; state<3; state++){
					for(int model=0; model<c.models; model++){
						scores[model]=state==0 ? legacyRanks[0].rawScores[model]
							: (state==1 ? queryRanks[0] : fixedRanks[0]).score(NcrnaDdlRanks.ALL_KEYS, model);
					}
					order[state].build(window.nativeCounts, scores);
				}
			}
			for(final Edge edge:window.edges){
				for(int model=0; model<c.models; model++){acceptable[model]=(edge.uniqueModels[model>>>6]&(1L<<(model&63)))!=0;}
				for(int method=0; method<TOTAL; method++){
					for(int state=0; state<3; state++){
						final int[] tuple=tuples[state];
						if(edge.unique==0){Arrays.fill(tuple, -1);}
						else if(method==0){order[state].nativeBest(acceptable, tuple);}
						else if(method==1){order[state].sketchBest(acceptable, tuple);}
						else if(method==2){order[state].hybridBest(acceptable, tuple);}
						else{
							final int control=method-3;
							final int g=1+control/2;
							final int score=(control&1)==0 ? NcrnaDdlRanks.SSU : NcrnaDdlRanks.ALL_KEYS;
							if(state==0){legacyRanks[g].best(score, acceptable, tuple);}
							else{(state==1 ? queryRanks[g] : fixedRanks[g]).best(score, acceptable, tuple);}
						}
						final int offset=(edge.parent.index*TOTAL+method)*3+state;
						minima[offset]=NcrnaDdlSavedData.minimum(minima[offset], tuple[0]);
					}
					if(method==0){
						require(Arrays.equals(tuples[0], edge.nativeTuple) && Arrays.equals(tuples[0], tuples[1]) && Arrays.equals(tuples[0], tuples[2]),
							"Original native tuples must reproduce exactly in every hash comparison state");
					}
					if(method>=3){
						for(int field=0; field<4; field++){require(tuples[0][field]==edge.controls[(method-3)*4+field], "Original canonical cached tuple must reproduce exactly");}
					}
					require(ns>0 || Arrays.equals(tuples[0], tuples[1]), "N-free queries cannot change against unchanged member profiles");
					final boolean queryChanged=!Arrays.equals(tuples[0], tuples[1]);
					final boolean modelChanged=!Arrays.equals(tuples[1], tuples[2]);
					if(tuples[0][0]!=tuples[2][0]){changedRanks[ns==0 ? 0 : 1][method]++;}
					if(queryChanged || modelChanged){
						changedTuples[method]++;
						row.clear().append(window.input.panel.nameBytes).tab().append(window.input.slot).tab().append(window.input.idBytes)
							.tab().append(window.event).tab().append(edge.parent.idBytes).tab().append(NcrnaDdlSavedData.STRATA[edge.parent.stratum])
							.tab().append(ns).tab().append(methodName(method));
						for(final int[] tuple:tuples){for(final int field:tuple){row.tab().append(field);}}
						batch.append(row.nl());
						if(batch.length()>=65536){flush();}
					}
				}
			}
		}
		private void flush(){
			if(batch.length()==0){return;}
			synchronized(changes.writer){changes.writer.print(batch);}
			batch.clear();
		}
		final int[] minima=new int[data.parents.size()*TOTAL*3];
		final long[] changedTuples=new long[TOTAL];
		final long[][] changedRanks=new long[2][TOTAL];
		final long[] windowCounts=new long[2];
		final LegacyRanks[] legacyRanks=new LegacyRanks[9];
		final NcrnaDdlRanks[] queryRanks=new NcrnaDdlRanks[9];
		final NcrnaDdlRanks[] fixedRanks=new NcrnaDdlRanks[9];
		final NcrnaDdlInterleave[] order=new NcrnaDdlInterleave[3];
		final float[] scores=new float[c.models];
		final boolean[] acceptable=new boolean[c.models];
		final int[][] tuples=new int[3][4];
		final ByteBuilder row=new ByteBuilder();
		final ByteBuilder batch=new ByteBuilder();
		final Output changes;
		boolean success;
	}

	private void report(final ArrayList<Worker> workers){
		final int[] minima=new int[data.parents.size()*TOTAL*3];
		Arrays.fill(minima, -1);
		final long[] tupleChanges=new long[TOTAL];
		final long[][] rankChanges=new long[2][TOTAL];
		final long[] counts=new long[2];
		for(final Worker worker:workers){
			require(worker.success, "Every replay worker must finish before reporting");
			for(int i=0; i<minima.length; i++){minima[i]=NcrnaDdlSavedData.minimum(minima[i], worker.minima[i]);}
			for(int m=0; m<TOTAL; m++){tupleChanges[m]+=worker.changedTuples[m];}
			for(int n=0; n<2; n++){
				for(int m=0; m<TOTAL; m++){rankChanges[n][m]+=worker.changedRanks[n][m];}
			}
			for(int n=0; n<2; n++){counts[n]+=worker.windowCounts[n];}
		}
		require(counts[0]+counts[1]==windows.length, "Every selected window, including unrankable ones, must be counted");
		final long[][] changed=new long[3][TOTAL];
		final long[][][] sums=new long[3][TOTAL][3];
		final long[] common=new long[3];
		try(Output loci=new Output(c.output+".loci.tsv.gz"); Output summary=new Output(c.output+".summary.tsv")){
			loci.writer.println("panel\tinput_slot\tinput_id\tparent\tseen_stratum\tmethod\told_rank\tquery_fixed_rank\tall_fixed_rank\tsource_disposition");
			final ByteBuilder row=new ByteBuilder();
			for(final Parent parent:data.parents){
				if(parent.nativeRank>0){common[parent.stratum]++;}
				for(int m=0; m<TOTAL; m++){
					final int offset=(parent.index*TOTAL+m)*3;
					final int expected=m<3 ? oldMin[parent.index*3+m] : parent.controlMin[m-3];
					require(minima[offset]==expected, "Every old saved locus minimum must reproduce, including unrankable cases");
					if(minima[offset]!=minima[offset+2]){changed[parent.stratum][m]++;}
					if(parent.nativeRank>0){
						for(int state=0; state<3; state++){
							require(minima[offset+state]>0, "Common loci must retain a complete ranking of acceptable models");
							sums[parent.stratum][m][state]+=minima[offset+state];
						}
					}
					row.clear().append(parent.input.panel.nameBytes).tab().append(parent.input.slot).tab().append(parent.input.idBytes)
						.tab().append(parent.idBytes).tab().append(NcrnaDdlSavedData.STRATA[parent.stratum]).tab().append(methodName(m));
					for(int state=0; state<3; state++){row.tab().append(minima[offset+state]);}
					loci.writer.print(row.tab().append(parent.reason).nl());
				}
			}
			summary.writer.println("seen_stratum\tmethod\tcommon_loci\told_rank_sum\tquery_fixed_rank_sum\tall_fixed_rank_sum\tchanged_locus_ranks_all_targets\tchanged_window_tuples_all_strata\tN_free_windows_all_strata\tN_containing_windows_all_strata\tchanged_N_free_window_ranks_all_strata\tchanged_N_containing_window_ranks_all_strata");
			for(int stratum=0; stratum<3; stratum++){
				for(int m=0; m<TOTAL; m++){
					row.clear().append(NcrnaDdlSavedData.STRATA[stratum]).tab().append(methodName(m)).tab().append(common[stratum]);
					for(int state=0; state<3; state++){row.tab().append(sums[stratum][m][state]);}
					summary.writer.print(row.tab().append(changed[stratum][m]).tab().append(tupleChanges[m]).tab().append(counts[0]).tab().append(counts[1])
						.tab().append(rankChanges[0][m]).tab().append(rankChanges[1][m]).nl());
				}
			}
		}
	}

	/** Preserved pre-fix loop, never called by the production caller. */
	private static void legacy(final byte[] bases, final int k, final boolean canonical, final boolean middle,
			final DynamicDemiLog sketch, final NcrnaDdlQueryKeys raw){
		final long mask=(1L<<(2*k))-1;
		final long middleMask=middle ? ~(3L<<(((k-1)/2)*2)) : -1L;
		final int lenmask=(-1)>>>k;
		long forward=0;
		long reverse=0;
		int len=-1;
		for(final byte b:bases){
			final long x=((b>>1)&3)|((b&8)<<28);
			len|=x;
			len>>>=1;
			forward=((forward<<2)|x)&mask;
			if(canonical){reverse=((reverse>>>2)|((x^2)<<(2*k-2)))&mask;}
			if(len>lenmask){continue;}
			final long word=(canonical ? Math.max(forward, reverse) : forward)&middleMask;
			if(sketch!=null){sketch.add(word);}
			if(raw!=null){raw.add(word);}
		}
	}

	private static final class LegacyRanks{
		LegacyRanks(final NcrnaDdlProfile[] profiles_, final int k_, final boolean canonical_, final boolean middle_){
			profiles=profiles_;
			k=k_;
			canonical=canonical_;
			middle=middle_;
			query=DynamicDemiLog.create(profiles[0].buckets(), k, 12345, 0, false, false);
			raw=new NcrnaDdlQueryKeys(profiles[0].buckets(), k);
			rawScores=new float[profiles.length];
			ssuScores=new float[profiles.length];
		}
		void score(final byte[] bases){
			query.clear();
			raw.build(EMPTY);
			legacy(bases, k, canonical, middle, query, raw);
			for(int m=0; m<profiles.length; m++){
				DynamicDemiLog.compareDetailed(query.maxArray(), profiles[m].keys(), counts);
				ssuScores[m]=DynamicDemiLog.wkid(counts[0], counts[1], counts[2])*(float)Math.sqrt(counts[1]);
				rawScores[m]=raw.equalBuckets(profiles[m].keys());
			}
		}
		void best(final int method, final boolean[] acceptable, final int[] tuple){
			NcrnaDdlRanks.best(method==NcrnaDdlRanks.SSU ? ssuScores : rawScores, acceptable, tuple);
		}
		final NcrnaDdlProfile[] profiles;
		final int k;
		final boolean canonical, middle;
		final DynamicDemiLog query;
		final NcrnaDdlQueryKeys raw;
		final float[] rawScores, ssuScores;
		final int[] counts=new int[4];
	}

	private static int countN(final byte[] bases){
		int count=0;
		for(final byte b:bases){
			if(b=='N' || b=='n'){count++;}
			else if(b!='A' && b!='C' && b!='G' && b!='T' && b!='a' && b!='c' && b!='g' && b!='t'){
				throw new IllegalArgumentException("The N-boundary replay requires ACGTN, not unsupported IUPAC byte="+b);
			}
		}
		return count;
	}
	private static int k(final int geometry){return geometry==0 ? 10 : geometry<=4 ? 9 : 11;}
	private static int width(final int geometry){return geometry==0 ? 2048 : WIDTHS[(geometry-1)&3];}
	private static String methodName(final int method){return method==0 ? "native" : method==1 ? "masked10_raw" : method==2 ? "hybrid" : "canonical_control_"+NcrnaDdlSavedData.CONTROL_METHODS[method-3];}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final class Config{
		public String panels, selection, baseline, members, ledger, assignments, census, summary, references, output;
		public int models=508, groups=3153, threads=16;
	}
	private final Config c;
	private NcrnaDdlSavedData data;
	private Window[] windows;
	private int[] oldMin;
	private final NcrnaDdlProfile[][] old=new NcrnaDdlProfile[9][];
	private final NcrnaDdlProfile[][] fresh=new NcrnaDdlProfile[9][];
	private final AtomicInteger next=new AtomicInteger();
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	private static final int TOTAL=19;
	private static final int[] WIDTHS={256,512,1024,2048};
	private static final byte[] EMPTY=new byte[0];
}
