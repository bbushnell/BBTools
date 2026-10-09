package assemble;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import map.LongHashMap;
import parse.Parse;
import parse.Parser;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Batch exact-reference labels, run only after read-count tables are released.
 * A query-seed index scans each reference orientation once and full-verifies hits.
 * Hit counts saturate at two: only absent, unique and multiple matter for labels.
 * @author Fischl
 */
final class FusionJoinLabeler {

	/** Parses native options after the launcher's PreParser has expanded its arguments. */
	static void run(final String[] args){
		final Parser parser=new Parser();
		String trace=null, reference=null;
		boolean verify=false;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase(java.util.Locale.ROOT);
			while(a.startsWith("-")){a=a.substring(1);}
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("trace")){trace=b;}
			else if(a.equals("ref")){reference=b;}
			else if(a.equals("verify")){verify=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown label option: "+arg);}
		}
		if(trace==null || reference==null || parser.out1==null){
			throw new IllegalArgumentException("Label mode requires trace= ref= out=; verify=t is a small-fixture exhaustive cross-check.");
		}
		if(!Tools.testInputFiles(false, true, trace, reference) ||
				!Tools.testOutputFiles(false, false, false, parser.out1)){
			throw new IllegalArgumentException("Label inputs must exist and output must be a new writable file.");
		}
		TadpoleGraph.checkPaths(reference, new ArrayList<String>(Arrays.asList(trace)), parser.out1);
		final ArrayList<FusionJoinDiagnostic.Join> joins=FusionJoinDiagnostic.readJoins(trace, true);
		final ArrayList<String> refs=FusionJoinDiagnostic.sequences(reference);
		final ArrayList<byte[]> patterns=new ArrayList<byte[]>(joins.size()*3);
		for(FusionJoinDiagnostic.Join join : joins){
			patterns.add(join.product);
			patterns.add(join.left);
			patterns.add(join.right);
		}
		final FusionJoinDiagnostic.Hit[] hits=findAll(refs, patterns);
		if(verify){verify(refs, patterns, hits);}
		writeLabels(joins, refs, hits, parser.out1);
	}

	/** Writes one truth record per proposal, including explicit unresolved cases. */
	private static void writeLabels(final ArrayList<FusionJoinDiagnostic.Join> joins,
			final ArrayList<String> refs, final FusionJoinDiagnostic.Hit[] hits, final String output){
		assert(hits.length==joins.size()*3) : "Each proposal requires whole, left and right match censuses.";
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
		writer.start();
		int supported=0, contradicted=0, unresolved=0;
		boolean error;
		try{
			writer.print("#counts_saturated_at=2; target=NA is excluded from binary training\n");
			writer.print("event\tphase_k\tquery_k\tlabel\ttarget\twhole_hits\tleft_hits\tright_hits\n");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<joins.size(); i++){
				final FusionJoinDiagnostic.Join join=joins.get(i);
				final FusionJoinDiagnostic.Hit whole=hits[3*i], left=hits[3*i+1], right=hits[3*i+2];
				final String label=FusionJoinDiagnostic.label(join, refs, whole, left, right);
				final String target;
				if(label.equals("supported")){target="1"; supported++;}
				else if(label.equals("contradicted")){target="0"; contradicted++;}
				else{target="NA"; unresolved++;}
				row.clear().append(join.event).tab().append(join.phase).tab().append(join.phase).tab().append(label);
				row.tab().append(target).tab().append(whole.count).tab().append(left.count).tab().append(right.count).nl();
				writer.print(row);
			}
		}finally{error=writer.poisonAndWait();}
		if(error){throw new RuntimeException("Fusion labels could not be written: "+output);}
		System.err.println("FUSION_LABEL_PASS rows="+joins.size()+" supported="+supported+
				" contradicted="+contradicted+" unresolved="+unresolved);
	}

	/**
	 * Indexes one diverse exact31mer per query. Saturated queries are unlinked so
	 * repeated reference words do not repeatedly verify already ambiguous patterns.
	 * Short/undefined-only queries use the exhaustive fallback; corpus windows are longer.
	 */
	static FusionJoinDiagnostic.Hit[] findAll(final ArrayList<String> refs, final ArrayList<byte[]> patterns){
		assert(refs!=null && patterns!=null) : "Batch truth lookup requires explicit reference records and patterns.";
		final LongHashMap index=new LongHashMap(Math.max(16, patterns.size()*2));
		final Query[] queries=new Query[patterns.size()];
		// LongHashMap.put is insert-if-absent, not replacement. Keep mutable list
		// heads in a primitive array and map each seed to an immutable group ID.
		final int[] heads=new int[patterns.size()];
		Arrays.fill(heads, -1);
		final FusionJoinDiagnostic.Hit[] hits=new FusionJoinDiagnostic.Hit[patterns.size()];
		if(patterns.isEmpty()){return hits;}
		for(int i=0; i<queries.length; i++){
			final Query q=new Query(patterns.get(i));
			queries[i]=q;
			hits[i]=q.hit;
			if(q.offset<0){
				final FusionJoinDiagnostic.Hit exact=FusionJoinDiagnostic.find(refs, q.bases);
				q.hit.count=Math.min(2, exact.count);
				q.hit.sequence=exact.sequence;
				q.hit.pos=exact.pos;
			}else{
				int group=index.get(q.seed);
				if(group<0){index.put(q.seed, i); group=i;}
				q.next=heads[group];
				heads[group]=i;
			}
		}
		for(int r=0; r<refs.size(); r++){
			final String ref=refs.get(r);
			long seed=0;
			int valid=0;
			for(int p=0; p<ref.length(); p++){
				final int base=AminoAcid.baseToNumber[ref.charAt(p)];
				if(base<0){valid=0; seed=0; continue;}
				seed=((seed<<2)|base)&MASK;
				if(++valid<SEED_LENGTH){continue;}
				final int group=index.get(seed);
				int previous=-1, id=group<0 ? -1 : heads[group];
				while(id>=0){
					final Query q=queries[id];
					final int next=q.next, start=p-SEED_LENGTH+1-q.offset;
					if(matches(ref, start, q.bases)){
						q.hit.count++;
						q.hit.sequence=r;
						q.hit.pos=start;
					}
					if(q.hit.count>=2){
						if(previous<0){heads[group]=next;}else{queries[previous].next=next;}
					}else{previous=id;}
					id=next;
				}
			}
		}
		return hits;
	}

	/** Seed equality only proposes a start; every query base must match inside one reference record. */
	private static boolean matches(final String ref, final int start, final byte[] query){
		assert(query.length>0) : "Empty truth queries would match every reference coordinate.";
		if(start<0 || (long)start+query.length>ref.length()){return false;}
		for(int i=0; i<query.length; i++){if(ref.charAt(start+i)!=(char)query[i]){return false;}}
		return true;
	}

	/** Independently compares indexed hit classes and unique coordinates with exhaustive substring lookup. */
	private static void verify(final ArrayList<String> refs, final ArrayList<byte[]> patterns,
			final FusionJoinDiagnostic.Hit[] hits){
		assert(patterns.size()==hits.length) : "Verification needs one result per original query.";
		for(int i=0; i<hits.length; i++){
			final FusionJoinDiagnostic.Hit expected=FusionJoinDiagnostic.find(refs, patterns.get(i)), actual=hits[i];
			if(actual.count!=Math.min(2, expected.count) || (actual.count==1 &&
					(actual.sequence!=expected.sequence || actual.pos!=expected.pos))){
				throw new AssertionError("Indexed truth mismatch at query "+i+": "+actual.count+" vs "+expected.count);
			}
		}
		System.err.println("FUSION_LABEL_EXHAUSTIVE_PARITY_PASS queries="+hits.length);
	}

	/** Shared-seed collisions, repeats, absent queries, short fallback and empty batches. */
	static void selfTest(){
		final String seed="ACGTACGTACGTACGTACGTACGTACGTACG";
		final ArrayList<byte[]> patterns=new ArrayList<byte[]>();
		for(String suffix : new String[]{"AAAAAAAAAA", "TTTTTTTTTT", "CCCCCCCCCC"}){
			patterns.add((seed+suffix).getBytes(StandardCharsets.US_ASCII));
		}
		patterns.add("GATTACA".getBytes(StandardCharsets.US_ASCII));
		final ArrayList<String> refs=new ArrayList<String>();
		refs.add(seed+"AAAAAAAAAANN"+seed+"TTTTTTTTTTNN"+seed+"AAAAAAAAAA");
		refs.add("GATTACANN");
		final FusionJoinDiagnostic.Hit[] hits=findAll(refs, patterns);
		verify(refs, patterns, hits);
		if(hits[0].count!=2 || hits[1].count!=1 || hits[2].count!=0 || hits[3].count!=1){
			throw new AssertionError("Shared-seed fixture did not distinguish full-query matches.");
		}
		if(findAll(refs, new ArrayList<byte[]>()).length!=0){throw new AssertionError("Empty batch produced truth records.");}
		System.err.println("FUSION_LABEL_TEST_PASS");
	}

	/** One persistent query plus primitive linked-list state; no per-reference-position allocation. */
	private static final class Query {
		/** Chooses the earliest31mer with the most balanced A/C/G/T composition. */
		Query(final byte[] bases_){
			assert(bases_!=null && bases_.length>0) : "Truth queries must contain sequence.";
			bases=bases_;
			final int[] counts=new int[4];
			int valid=0, best=Integer.MAX_VALUE;
			long word=0;
			for(int p=0; p<bases.length; p++){
				final int base=AminoAcid.baseToNumber[bases[p]];
				if(base<0){valid=0; word=0; Arrays.fill(counts, 0); continue;}
				word=((word<<2)|base)&MASK;
				counts[base]++;
				if(++valid>SEED_LENGTH){counts[AminoAcid.baseToNumber[bases[p-SEED_LENGTH]]]--;}
				if(valid<SEED_LENGTH){continue;}
				int score=0;
				for(int n : counts){score+=n*n;}
				if(score<best){best=score; offset=p-SEED_LENGTH+1; seed=word;}
			}
		}
		final byte[] bases;
		final FusionJoinDiagnostic.Hit hit=new FusionJoinDiagnostic.Hit();
		long seed;
		int offset=-1, next=-1;
	}

	private static final int SEED_LENGTH=31;
	private static final long MASK=(1L<<(2*SEED_LENGTH))-1;
}
