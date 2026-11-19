package prok;

import java.util.Arrays;

import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Optional exact-sequence partition of already validated native grade results.
 * The benchmark/subtype consumer owns matching and outcome validation. This
 * helper joins original escaped input/parent bytes, never coordinates, and
 * preserves UNKNOWN separately from NOT_SEEN. Query labels describe their
 * matched reference, not an independently measured query sequence. Unmatched
 * and out-of-population queries remain UNASSIGNED and precision remains NA.
 * Native byte I/O follows template.A_SampleByteFile; only retained keys allocate.
 * @author Ganyu
 */
final class NcrnaExactMembership{

	NcrnaExactMembership(final String path, final ObjectIntMap<Target> parents,
			final int parentCount, final byte[][] contexts_){
		require(path!=null && parents!=null && contexts_!=null && contexts_.length>0,
			"Exact membership needs the validated parent universe and declared cells");
		parentIndex=parents;
		contexts=contexts_;
		exact=new byte[parentCount];
		source=new byte[parentCount];
		complete=new boolean[parentCount];
		common=new boolean[parentCount];
		Arrays.fill(exact, (byte)-1);
		counts=new long[contexts.length*3*3*3];
		for(int cell=0; cell<contexts.length; cell++){
			probe.set(contexts[cell], 0, contexts[cell].length);
			require(contextIndex.get(probe)<0, "Declared exact-membership cells must be unique");
			contextIndex.put(probe.copy(), cell);
		}
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaDdlCallerPanel.STRATA_HEADER, path);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==6, "Exact membership requires the six-field caller-panel strata schema");
				final int id=parent(p, 0, 1);
				require(id>=0, "Exact label has no original input/parent in the validated union grades");
				require(exact[id]<0, "Duplicate exact label for an original input/parent");
				exact[id]=(byte)state(2);
				source[id]=(byte)state(3);
				complete[id]=bool(4);
				common[id]=bool(5);
				require(complete[id]==(exact[id]!=2), "Incomplete/undefined sequence must be UNKNOWN; complete sequence must have a definite exact label");
				labels++;
			}
		}
	}

	boolean labelled(final int id){return id>=0 && id<exact.length && exact[id]>=0;}
	int size(){return labels;}

	/** outcome uses the native consumer's TP=1,FN=2,DONT_CARE=3 contract. */
	void reference(final int cell, final int view, final int id, final byte[] key,
			final byte[] geometry, final int outcome, final long queryLine, final Output out){
		require(labelled(id) && cell>=0 && cell<contexts.length && view>=0 && view<3
			&& outcome>=1 && outcome<=3, "Only bound labelled native reference outcomes enter exact-state counts");
		require((outcome==1)==(queryLine>0), "A labelled TP must retain its original matched query-line witness");
		final int offset=((cell*3+view)*3+exact[id])*3+outcome-1;
		counts[offset]=Math.addExact(counts[offset], 1);
		b.clear().append(contexts[cell]).tab().append(VIEWS[view]).tab().append(key).tab()
			.append(STATES[exact[id]]).tab().append(STATES[source[id]]).tab().append(complete[id]).tab().append(common[id])
			.tab().append(geometry).tab().append(OUTCOMES[outcome-1]).tab();
		if(queryLine>0){b.append(queryLine);}
		else{b.append('.');}
		out.writer.print(b.nl());
	}

	/** Re-read only validated text to retain every original query row verbatim.
	 * No per-query row/string is retained in heap. Counts must still reconcile
	 * with the complete query populations held by the validating consumer. */
	void finish(final String prefix, final String unionRaw, final String subtypeRaw,
			final long[] expectedQueries){
		require(expectedQueries.length==contexts.length*3, "Every cell/view needs its validated query population");
		final long[][] queries=new long[expectedQueries.length][7];
		queryRows(unionRaw, prefix+SUFFIXES[2], false, queries);
		queryRows(subtypeRaw, prefix+SUFFIXES[3], true, queries);
		try(Output out=new Output(prefix+SUFFIXES[4])){
			out.writer.println("policy\tthreads\trepeat\tview\tquery_rows\tTP\tFP\tMASKED\tEXCLUDED\tduplicate_queries\tunassigned_reference_state\tqualification");
			for(int cell=0; cell<contexts.length; cell++){
				for(int view=0; view<3; view++){
					final int index=cell*3+view;
					final long[] q=queries[index];
					require(q[0]==expectedQueries[index] && q[0]==q[1]+q[2]+q[3]+q[4],
						"Exact query witnesses must conserve every validated original query, including unmatched and duplicate rows");
					b.clear().append(contexts[cell]).tab().append(VIEWS[view]);
					for(final long value:q){b.tab().append(value);}
					out.writer.print(b.append("\tQUERY_ROWS_NOT_UNIQUE_LOCI;UNION_INCLUDES_ALL_COBALAMIN_MODELS\n"));
				}
			}
		}
		try(Output out=new Output(prefix+SUFFIXES[0])){
			out.writer.println("policy\tthreads\trepeat\tview\texact_sequence\tTP\tFN\tmasked_reference_parents\tscored_reference_parents\tprecision\trecall\tqualification");
			for(int cell=0; cell<contexts.length; cell++){
				for(int view=0; view<3; view++){
					for(int state=0; state<=3; state++){
						long tp=0;
						long fn=0;
						long masks=0;
						for(int s=0; s<3; s++){
							if(state<3 && s!=state){continue;}
							final int offset=((cell*3+view)*3+s)*3;
							tp=Math.addExact(tp, counts[offset]);
							fn=Math.addExact(fn, counts[offset+1]);
							masks=Math.addExact(masks, counts[offset+2]);
						}
						final long scored=Math.addExact(tp, fn);
						b.clear().append(contexts[cell]).tab().append(VIEWS[view]).tab().append(state==3 ? "ALL_LABELLED" : STATES[state])
							.tab().append(tp).tab().append(fn).tab().append(masks).tab().append(scored).append("\tNA\t");
						if(scored==0){b.append("NA");}
						else{b.append(tp/(double)scored, 6);}
						out.writer.print(b.append("\tRF00174_LABELLED_PARENTS;PRECISION_UNAVAILABLE_UNMATCHED_QUERIES_UNPARTITIONED\n"));
					}
				}
			}
		}
	}

	private void queryRows(final String path, final String output, final boolean subtype,
			final long[][] queries){
		final String header="policy\tthreads\trepeat\t"+(subtype ? NcrnaSubtypeGrade.RAW_HEADER : GradeNcrnaReference.RAW_HEADER);
		try(Input in=new Input(path); Output out=new Output(output)){
			RrnaResourceIO.header(in.reader, header, path);
			out.writer.println("reference_exact_sequence\t"+header);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==(subtype ? 25 : 21), "Native query witness schema changed after validation");
				if(!p.termEquals("QUERY", subtype ? 6 : 4)){continue;}
				if(!p.termEquals(subtype ? "RF00174" : "cobalamin", subtype ? 4 : 5)){continue;}
				int view=0;
				if(subtype){
					view=p.termEquals(NcrnaSubtypeGrade.VIEWS[0], 3) ? 1 : p.termEquals(NcrnaSubtypeGrade.VIEWS[1], 3) ? 2 : -1;
					require(view>0, "Unknown native subtype view cannot enter exact query witnesses");
				}
				key.clear();
				for(int i=0; i<3; i++){
					if(i>0){key.tab();}
					field(key, p, i);
				}
				probe.set(key.array, 0, key.length());
				final int cell=contextIndex.get(probe);
				require(cell>=0, "Query witness has an undeclared benchmark cell");
				final long[] q=queries[cell*3+view];
				q[0]++;
				final int status=subtype ? 14 : 10;
				int label=-1;
				if(p.termEquals("TP", status)){
					q[1]++;
					final int id=parent(p, subtype ? 5 : 3, subtype ? 15 : 13);
					require(id>=0, "A validated TP must retain its original matched parent");
					if(labelled(id)){label=exact[id];}
				}else if(p.termEquals("FP", status)){q[2]++;}
				else if(p.termEquals("MASKED", status)){q[3]++;}
				else{
					require(!subtype && p.termEquals("EXCLUDED", status), "Unknown query outcome cannot disappear from the saved population");
					q[4]++;
				}
				if(p.termEquals("true", subtype ? 18 : 16)){q[5]++;}
				if(label<0){q[6]++;}
				out.writer.print(b.clear().append(label<0 ? "UNASSIGNED" : STATES[label]).tab().append(row).nl());
			}
		}
	}

	private int parent(final LineParser1 parser, final int input, final int parent){
		key.clear();
		field(key, parser, input).tab();
		field(key, parser, parent);
		probe.set(key.array, 0, key.length());
		return parentIndex.get(probe);
	}
	private int state(final int field){
		for(int i=0; i<STATES.length; i++){if(p.termEquals(STATES[i], field)){return i;}}
		throw new IllegalArgumentException("Unknown exact-sequence/source-locus state; UNKNOWN must stay explicit");
	}
	private boolean bool(final int field){
		require(p.termEquals("true", field) || p.termEquals("false", field), "Membership booleans must be literal true/false");
		return p.termEquals("true", field);
	}
	private static ByteBuilder field(final ByteBuilder b, final LineParser1 p, final int field){
		p.setBounds(field);
		return b.append(p.line(), p.a(), p.b()-p.a());
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	private final ObjectIntMap<Target> parentIndex;
	private final ObjectIntMap<Target> contextIndex=new ObjectIntMap<Target>(Target.class);
	private final byte[][] contexts;
	private final byte[] exact, source;
	private final boolean[] complete, common;
	private final long[] counts;
	private int labels;
	private final LineParser1 p=new LineParser1('\t');
	private final Target probe=new Target();
	private final ByteBuilder key=new ByteBuilder();
	private final ByteBuilder b=new ByteBuilder();
	static final String[] SUFFIXES={".exact_summary.tsv", ".exact_references.tsv", ".exact_union_queries.tsv", ".exact_subtype_queries.tsv", ".exact_query_counts.tsv"};
	private static final String[] STATES={"SEEN", "NOT_SEEN", "UNKNOWN"};
	private static final String[] OUTCOMES={"TP", "FN", "DONT_CARE"};
	private static final String[] VIEWS={"COBALAMIN_UNION_ON_RF00174_PARENTS", "MODEL_SUPPORT_RECOVERY", "EXCLUSIVE_SUBTYPE_ACCURACY"};
	static final String REFERENCE_HEADER="policy\tthreads\trepeat\tview\tinput_id\treference_id\texact_sequence\tsame_source_locus\tcomplete_defined_query\tcommon_native_rankable\tseqid\tstart\tstop\tstrand\toutcome\tmatched_query_line";
}
