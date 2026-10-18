package test.prot;

import prot.DigestSuffix;
import prot.HbmCompetitiveAssign;
import prot.HmmComparisonData;
import prot.ProteinFamilyAssigner;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/** Compares the reusable API with frozen experimental calls, including threading.
 * Runs on the cluster with the complete reference; not part of laptop fixtures.
 * @author Keqing
 */
public final class ProteinFamilyAssignerTest {

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, ProteinFamilyAssignerTest.class, false);
		try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final HashMap<String,String> options) throws Exception{
		for(String key : options.keySet()){
			require(Arrays.asList("models", "in", "expected", "count", "t", "out").contains(key), "Unknown API test parameter: "+key);
		}
		final int count=Integer.parseInt(HmmComparisonData.required(options, "count"));
		final int threads=Integer.parseInt(options.getOrDefault("t", "4"));
		require(count>0 && threads>1 && threads<=32, "API parity test requires a positive count and 2..32 threads");
		final Path output=Paths.get(HmmComparisonData.required(options, "out"));
		require(!Files.exists(output), "API test output must be a fresh directory");
		Shared.setThreads(threads); Shared.AMINO_IN=true;
		final ProteinFamilyAssigner library=ProteinFamilyAssigner.load(Paths.get(HmmComparisonData.required(options, "models")));
		final ArrayList<Query> queries=readQueries(HmmComparisonData.required(options, "in"));
		require(queries.size()==count, "Frozen query count differs: "+queries.size()+" != "+count);
		final ArrayList<byte[]> expected=ByteFile.toLines(HmmComparisonData.required(options, "expected"));
		require(expected.size()==count+1 && HEADER.equals(new String(expected.get(0), java.nio.charset.StandardCharsets.US_ASCII)), "Frozen assignment header/count differs");
		final ProteinFamilyAssigner.Hit[] single=new ProteinFamilyAssigner.Hit[count], parallel=new ProteinFamilyAssigner.Hit[count];
		final ProteinFamilyAssigner.Worker worker=library.newWorker();
		for(int i=0; i<count; i++){final Query q=queries.get(i); single[i]=worker.assign(q.id, q.raw);}
		final AtomicInteger next=new AtomicInteger(); final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final ArrayList<Thread> workers=new ArrayList<Thread>();
		for(int t=0; t<threads; t++){
			workers.add(new Thread(()->{
				try{
					final ProteinFamilyAssigner.Worker local=library.newWorker();
					for(int start=next.getAndAdd(8); start<count && failure.get()==null; start=next.getAndAdd(8)){
						for(int i=start, stop=Math.min(count, start+8); i<stop; i++){
							final Query q=queries.get(i); parallel[i]=local.assign(q.id, q.raw);
						}
					}
				}catch(Throwable e){failure.compareAndSet(null, e);}
			}));
		}
		ThreadWaiter.startAndWait(workers);
		if(failure.get()!=null){throw new RuntimeException("Threaded family API test failed", failure.get());}
		Files.createDirectory(output);
		final ByteStreamWriter one=HmmComparisonData.writer(output.resolve("single.tsv").toString());
		final ByteStreamWriter many=HmmComparisonData.writer(output.resolve("parallel.tsv").toString());
		final ByteBuilder row=new ByteBuilder(); long assigned=0;
		try{
			one.println(HEADER); many.println(HEADER);
			for(int i=0; i<count; i++){
				final byte[] a=format(row, queries.get(i).id, single[i]), b=format(row, queries.get(i).id, parallel[i]);
				require(Arrays.equals(a, expected.get(i+1)), "API differs from frozen call for query "+queries.get(i).id);
				require(Arrays.equals(a, b), "Thread count changed the call for query "+queries.get(i).id);
				one.println(a); many.println(b); if(single[i].assigned()){assigned++;}
			}
		}finally{HmmComparisonData.close(one); HmmComparisonData.close(many);}
		final ByteStreamWriter receipt=HmmComparisonData.writer(output.resolve("PASS").toString());
		receipt.println("FAMILY_API_PARITY_PASS queries="+count+" assigned="+assigned+" threads="+threads+" frozen_exact=true threaded_exact=true");
		receipt.println("model_sha80="+library.modelPin());
		receipt.println("identity_sha80="+library.identityPin());
		receipt.println("consensus_sha80="+library.consensusPin());
		receipt.println("queries_sha80="+DigestSuffix.file(options.get("in")));
		receipt.println("expected_sha80="+DigestSuffix.file(options.get("expected")));
		HmmComparisonData.close(receipt);
		System.err.println("FAMILY_API_PARITY_PASS queries="+count+" assigned="+assigned+" threads="+threads);
	}

	private static ArrayList<Query> readQueries(final String path){
		final ArrayList<Query> queries=new ArrayList<Query>();
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA, null, true, false), null, true, -1);
		input.start();
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null && batch.size()>0; batch=input.nextList()){
				for(Read read : batch){
					require(read.mate==null, "Frozen protein queries must be unpaired");
					final String id=HbmCompetitiveAssign.queryId(read.id);
					//Match the frozen experimental CLI's explicit malformed-edge repair.
					queries.add(new Query(id, HbmCompetitiveAssign.trimEdges(read.bases, id)));
				}
			}
		}finally{input.close();}
		require(!input.errorState(), "Frozen protein query read failed");
		return queries;
	}

	private static byte[] format(final ByteBuilder row, final String id, final ProteinFamilyAssigner.Hit hit){
		require(hit!=null, "A worker failed to produce a query result"); row.clear().append(id).tab();
		if(hit.assigned()){
			row.append("ASSIGNED").tab().append(hit.rank).tab().append(hit.familyId).tab().append(hit.representative).tab().append(hit.score64)
				.tab().append(hit.identity, 8).tab().append(hit.queryCoverage, 8).tab().append(hit.coreCoverage, 8).tab().append(hit.mutualCoverage, 8)
				.tab().append(hit.queryLength).tab().append(hit.referenceStart).tab().append(hit.referenceEnd);
		}else{row.append("NO_ELIGIBLE_FAMILY\t-1\t-1\tNA\tNA\tNA\tNA\tNA\tNA").tab().append(hit.queryLength).append("\t-1\t-1");}
		return row.toBytes();
	}
	private static void require(boolean ok, String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	private static final class Query {
		Query(String id_, byte[] raw_){id=id_; raw=raw_;}
		final String id;
		final byte[] raw;
	}
	private static final String HEADER="query\tstatus\tactive_index\tfamily_id\trep_id\tscore64\tidentity_pct\tcoverage_q\tcoverage_core\tmutual_coverage\tquery_length\tref_start\tref_end";
}
