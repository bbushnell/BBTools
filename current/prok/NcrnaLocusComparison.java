package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import structures.IntList;

/** One-to-one detection matching for a declared ncRNA family.
 * Coordinates are one-based inclusive. Compatible calls overlap at least half
 * the shorter span on the same contig and strand. Augmenting paths find a
 * maximum-cardinality matching; tied assignments are not a biological claim.
 * Each locus is a retained result object. Primitive arrays hold matching state.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLocusComparison{

	/** Reads the five-column truth table, retaining only the exact requested label. */
	public static ArrayList<Locus> readTruth(final String path, final String label){
		final ArrayList<Locus> result=new ArrayList<Locus>();
		final LineParser1 lp=new LineParser1('\t');
		try(Input input=new Input(path)){
			final ByteFile file=input.reader;
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lp.terms()!=5){throw new IllegalArgumentException("Truth requires seqid, family, start, stop, strand: "+path);}
				if(lp.termEquals("seqid", 0) && lp.termEquals("family", 1)){continue;}
				if(!lp.termEquals(label, 1)){continue;}
				result.add(new Locus(lp.parseString(0), position(lp, 2), position(lp, 3), strand(lp, 4), null));
			}
		}
		if(result.isEmpty()){throw new IllegalArgumentException("No reference loci have truthfamily="+label);}
		return result;
	}

	/** Reads ordinary caller GFF and retains the exact biological family label.
	 * The caller's model attribute supplies original-winner restoration provenance.
	 */
	public static ArrayList<Locus> readCalls(final String path, final String label){
		final ArrayList<Locus> result=new ArrayList<Locus>();
		final LineParser1 lp=new LineParser1('\t');
		try(Input input=new Input(path)){
			final ByteFile file=input.reader;
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lp.terms()!=9){throw new IllegalArgumentException("Caller output must be nine-column GFF: "+path);}
				lp.setBounds(8);
				final int from=lp.a(), to=lp.b();
				final int familyStart=attributeStart(line, from, to, "family:");
				if(familyStart<0){continue;}
				final int familyEnd=attributeEnd(line, familyStart, to);
				if(!equals(line, familyStart, familyEnd, label)){continue;}
				final int modelStart=attributeStart(line, from, to, "model:");
				if(modelStart<0){throw new IllegalArgumentException("Selected caller family lacks a winning model attribute");}
				int modelEnd=attributeEnd(line, modelStart, to);
				for(int i=modelStart; i<modelEnd; i++){
					if(line[i]==' '){
						modelEnd=i;
						break;
					}
				}
				final String model=new String(line, modelStart, modelEnd-modelStart, StandardCharsets.US_ASCII);
				result.add(new Locus(lp.parseString(0), position(lp, 3), position(lp, 4), strand(lp, 6), model));
			}
		}
		return result;
	}

	/** Returns one call index per reference, or -1 if unmatched. Inputs are borrowed. */
	public static int[] match(final ArrayList<Locus> truth, final ArrayList<Locus> calls){
		return match(truth, calls, null);
	}

	/** Maximizes protected-reference matches first, then total matches.
	 * Augmenting the resulting matching never unmatches an already matched
	 * reference, so caller output order cannot invent an avoidable protected loss.
	 * A null priority array requests ordinary maximum-cardinality matching.
	 */
	public static int[] match(final ArrayList<Locus> truth, final ArrayList<Locus> calls,
			final boolean[] protectedReferences){
		assert(truth!=null && calls!=null) : "Detection matching requires explicit reference and prediction populations";
		if(protectedReferences!=null && protectedReferences.length!=truth.size()){
			throw new IllegalArgumentException("Protected flags must use the complete truth index axis");
		}
		final ObjectMap<String, IntList> byContig=new ObjectMap<String, IntList>(String.class, IntList.class);
		for(int reference=0; reference<truth.size(); reference++){
			final String name=truth.get(reference).seqid;
			IntList list=byContig.get(name);
			if(list==null){
				list=new IntList();
				byContig.put(name, list);
			}
			list.add(reference);
		}
		final IntList[] compatible=new IntList[calls.size()];
		for(int call=0; call<calls.size(); call++){
			final Locus prediction=calls.get(call);
			final IntList references=byContig.get(prediction.seqid);
			if(references==null){continue;}
			for(int i=0; i<references.size; i++){
				final int reference=references.get(i);
				if(overlaps(truth.get(reference), prediction)){
					if(compatible[call]==null){compatible[call]=new IntList();}
					compatible[call].add(reference);
				}
			}
		}
		final int[] referenceToCall=new int[truth.size()];
		Arrays.fill(referenceToCall, -1);
		final int[] visited=new int[truth.size()];
		// One reusable frame per possible call in an augmenting path. A long
		// overlap chain must not consume the JVM stack in proportion to its size.
		final int[] callStack=new int[calls.size()];
		final int[] nextOption=new int[calls.size()];
		final int[] viaReference=new int[calls.size()];
		final boolean[] used=protectedReferences==null ? null : new boolean[calls.size()];
		if(protectedReferences!=null){
			for(int call=0; call<calls.size(); call++){
				augment(call, compatible, referenceToCall, visited, call+1, protectedReferences,
					callStack, nextOption, viaReference);
			}
			for(final int call:referenceToCall){if(call>=0){used[call]=true;}}
			// The second pass reuses call-based visit IDs in a different search.
			Arrays.fill(visited, 0);
		}
		for(int call=0; call<calls.size(); call++){
			if(used==null || !used[call]){
				augment(call, compatible, referenceToCall, visited, call+1, null,
					callStack, nextOption, viaReference);
			}
		}
		return referenceToCall;
	}

	/** Depth-first augmentation using reused primitive frames instead of recursive calls.
	 * Each frame resumes at its next option after a failed child. On finding a free
	 * reference, unwind the path and replace its matches in the original visit order.
	 */
	private static boolean augment(final int call, final IntList[] compatible,
			final int[] assigned, final int[] visited, final int visit, final boolean[] eligible,
			final int[] callStack, final int[] nextOption, final int[] viaReference){
		assert(call>=0 && call<compatible.length) : "Augmenting path starts on the call axis: call="+call+", calls="+compatible.length;
		int depth=0;
		callStack[0]=call;
		nextOption[0]=0;
		while(depth>=0){
			final IntList options=compatible[callStack[depth]];
			if(options==null || nextOption[depth]>=options.size){
				depth--;
				continue;
			}
			final int reference=options.get(nextOption[depth]++);
			if(eligible!=null && !eligible[reference]){continue;}
			if(visited[reference]==visit){continue;}
			visited[reference]=visit;
			if(assigned[reference]<0){
				assigned[reference]=callStack[depth];
				while(depth>0){
					depth--;
					assigned[viaReference[depth]]=callStack[depth];
				}
				return true;
			}
			assert(depth+1<callStack.length)
				: "One-to-one matching and visited references prevent repeated calls on a path: depth="+depth+", calls="+callStack.length;
			viaReference[depth]=reference;
			depth++;
			callStack[depth]=assigned[reference];
			nextOption[depth]=0;
		}
		return false;
	}

	private static boolean overlaps(final Locus a, final Locus b){
		if(a.strand!=b.strand){return false;}
		final long overlap=(long)Math.min(a.stop, b.stop)-Math.max(a.start, b.start)+1;
		return overlap>0 && 2*overlap>=Math.min((long)a.stop-a.start+1, (long)b.stop-b.start+1);
	}

	/** Coordinates fit signed int after validation; intermediate overlap lengths use long. */
	private static int position(final LineParser1 lp, final int term){
		lp.setBounds(term);
		if(lp.b()-lp.a()<1 || lp.b()-lp.a()>10){throw new IllegalArgumentException("Invalid coordinate width");}
		final int position=NcrnaNumbers.unsignedInt(lp.line(), lp.a(), lp.b());
		if(position<1){throw new IllegalArgumentException("Coordinates must be positive: "+position);}
		return position;
	}

	private static int strand(final LineParser1 lp, final int term){
		if(lp.termEquals('+', term)){return 0;}
		if(lp.termEquals('-', term)){return 1;}
		throw new IllegalArgumentException("Detection comparison requires a known + or - strand");
	}

	private static int attributeStart(final byte[] line, final int from, final int to, final String key){
		for(int i=from; i+key.length()<=to; i++){
			if(i>from && line[i-1]!=',' && line[i-1]!=';'){continue;}
			if(equals(line, i, i+key.length(), key)){return i+key.length();}
		}
		return -1;
	}

	private static int attributeEnd(final byte[] line, final int from, final int to){
		int end=from;
		while(end<to && line[end]!=',' && line[end]!=';'){end++;}
		return end;
	}

	private static boolean equals(final byte[] line, final int from, final int to, final String value){
		if(to-from!=value.length()){return false;}
		for(int i=0; i<value.length(); i++){
			if(line[from+i]!=value.charAt(i)){return false;}
		}
		return true;
	}

	/** A retained reference or caller record; model is null for truth. */
	public static final class Locus{
		/** Retains an immutable one-based inclusive span on a named contig and known strand. */
		public Locus(final String seqid_, final int start_, final int stop_, final int strand_, final String model_){
			if(seqid_==null || seqid_.isEmpty() || start_<1 || stop_<start_ || (strand_!=0 && strand_!=1)){
				throw new IllegalArgumentException("Invalid locus coordinates, contig or strand");
			}
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			model=model_;
		}
		/** Exact contig identifier used to join calls and references. */
		public final String seqid;
		/** Winning caller model, or null for truth records. */
		public final String model;
		/** One-based inclusive genomic endpoints. */
		public final int start, stop;
		/** Zero for plus, one for minus. */
		public final int strand;
	}

	private NcrnaLocusComparison(){}
}
