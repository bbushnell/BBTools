package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import idaligner.Tracer;
import map.ObjectIntMap;
import map.ObjectSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

/** Compares every annotated mitochondrial rRNA to both universal SSU consensuses.
 * Exact-sequence caching preserves every original locus through a member map.
 * Shorter sequence is the query, as in NcrnaMitoSsuCorpus; no reverse-complement
 * search, score filter or biological classification. Coverage is nongap paired
 * columns divided by each full input length, not proof of an intact SSU.
 * ref_start0/ref_stop0 are zero-based inclusive coordinates in the longer input;
 * query_side identifies which input is the complete shorter query.
 * @author Brian Bushnell, Raiden */
public final class NcrnaMitoSsuCompare implements Accumulator<NcrnaMitoSsuCompare.Worker> {

	public static void main(String[] args){
		final NcrnaMitoSsuCompare tool=new NcrnaMitoSsuCompare(args);
		try{tool.process();}finally{Shared.closeStream(tool.outstream);}
	}
	private NcrnaMitoSsuCompare(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);outstream=pp.outstream;
		final Parser parser=new Parser();
		for(String arg:pp.args){
			final String[] parts=Parse.splitOnFirst(arg, '=');
			final String key=parts[0].toLowerCase(java.util.Locale.ROOT), value=parts.length>1 ? parts[1] : null;
			if(key.equals("ref16")){ref16=value;}
			else if(key.equals("ref18")){ref18=value;}
			else if(key.equals("outprefix")){prefix=value;}
			else if(key.equals("expected")){expected=Parse.parseIntKMG(value);}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		input=parser.in1;
		require(input!=null && ref16!=null && ref18!=null && prefix!=null && expected>0,
			"Require in= ref16= ref18= outprefix= expected=N>0");
		require(!parser.append && parser.maxReads<0, "Appending or limiting reads would break complete locus conservation");
		require(!Tracer.GLOBAL, "NcrnaMitoSsuCorpus comparison uses the default semiglobal Quantum/Tracer policy");
		Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;ReadWrite.setZipThreads(Math.min(4, Shared.threads()));
	}
	private void process(){
		for(String suffix:SUFFIXES){require(!new File(prefix+suffix).exists(), "Refusing existing comparison output: "+prefix+suffix);}
		refs[0]=reference(ref16, "universal");refs[1]=reference(ref18, "euk18S_universal");
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		stream(input, r->{
			final String name=first(r.id);require(names.add(name), "Duplicate source locus: "+name);
			final Key key=new Key(bases(r.bases));int index=cache.get(key);
			if(index<0){index=unique.size();cache.put(key, index);unique.add(key);}
			unique.get(index).copies++;members.add(new Member(name, index));
		});
		require(members.size()==expected, "Input locus count differs from expected: "+members.size()+" versus "+expected);
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		for(int i=0; i<Math.min(Shared.threads(), unique.size()); i++){workers.add(new Worker(this));}
		require(ThreadWaiter.startAndWait(workers, this) && !error, "Alignment worker failed; no comparison output is accepted");
		writeResults();
		outstream.println("MITO_SSU_COMPARE_PASS loci="+members.size()+" unique="+unique.size()+" alignments="+(2L*unique.size()));
	}
	private void writeResults(){
		final ByteStreamWriter rows=writer(prefix+SUFFIXES[0]), map=writer(prefix+SUFFIXES[1]);
		final Writer fasta=WriterFactory.makeWriter(FileFormat.testOutput(prefix+SUFFIXES[2], FileFormat.FASTA, null, true, false, false, true), null, 1);
		final ByteBuilder b=new ByteBuilder();fasta.start();
		final long[][] counts=new long[2][3];final double[][] sums=new double[2][2];long copies=0;
		try{
			rows.println("unique_id\tmodel\tcandidate_length\tmodel_length\tquery_side\tquantum_identity\ttrace_identity\tref_start0\tref_stop0\tmatches\tsubstitutions\tinsertions\tdeletions\tneutral_N\tpaired_columns\tcandidate_paired_fraction\tmodel_paired_fraction");
			map.println("locus_id\tunique_id\tlength");
			for(Member m:members){map.print(b.clear().append(m.name).tab().append("seq_").append(m.index).tab().append(unique.get(m.index).bases.length).nl());}
			ArrayList<Read> batch=new ArrayList<Read>(256);long batchId=0;
			for(int i=0; i<unique.size(); i++){
				final Key key=unique.get(i);copies+=key.copies;
				batch.add(new Read(key.bases, null, "seq_"+i, i, false));
				if(batch.size()==256){fasta.add(batch, batchId++);batch=new ArrayList<Read>(256);}
				for(int model=0; model<2; model++){
					final Result r=key.results[model];require(r!=null, "Every unique sequence must have both alignments");
					rows.print(b.clear().append("seq_").append(i).tab().append(MODELS[model]).tab().append(key.bases.length).tab().append(refs[model].length)
						.tab().append(r.candidateQuery ? "candidate" : "model").tab().append(r.identity, 9).tab().append(r.traceIdentity, 9)
						.tab().append(r.start).tab().append(r.stop).tab().append(r.matches).tab().append(r.subs).tab().append(r.ins).tab().append(r.dels)
						.tab().append(r.ns).tab().append(r.paired).tab().append((double)r.paired/key.bases.length, 9).tab().append((double)r.paired/refs[model].length, 9).nl());
					sums[0][model]+=key.copies*(double)r.identity;sums[1][model]+=r.identity;
				}
				final float a=key.results[0].identity, c=key.results[1].identity;
				final int winner=a>c ? 0 : a==c ? 1 : 2;counts[0][winner]+=key.copies;counts[1][winner]++;
			}
			if(!batch.isEmpty()){fasta.add(batch, batchId);}
		}finally{
			final boolean failed=rows.poisonAndWait()|map.poisonAndWait()|fasta.poisonAndWait();
			require(!failed, "Comparison/member/unique FASTA output failed");
		}
		require(copies==expected && fasta.readsWritten()==unique.size(), "Member and unique FASTA counts must conserve the full input population");
		final ByteStreamWriter summary=writer(prefix+SUFFIXES[3]);
		try{
			summary.println("weighting\tsequences\tcloser_16S\tequal\tcloser_18S\tmean_quantum_16S\tmean_quantum_18S\tunique_sequences\talignments");
			for(int w=0; w<2; w++){
				final int n=w==0 ? expected : unique.size();
				require(counts[w][0]+counts[w][1]+counts[w][2]==n, "Nearest-reference preferences must conserve each weighting denominator");
				summary.print(b.clear().append(w==0 ? "copy" : "unique").tab().append(n).tab().append(counts[w][0]).tab().append(counts[w][1])
					.tab().append(counts[w][2]).tab().append(sums[w][0]/n, 9).tab().append(sums[w][1]/n, 9).tab().append(unique.size()).tab().append(2L*unique.size()).nl());
			}
		}finally{require(!summary.poisonAndWait(), "Comparison summary write failed");}
	}
	/** QuantumAligner:535 enforces full-query and reference-span operation conservation.
	 * Returned identity retains neutral-N half credit; stats.identity is trace M/(M+S+I+D+N). */
	public static Result compare(byte[] candidate, byte[] model, AlignmentStats stats){
		require(candidate!=null && model!=null && candidate.length>0 && model.length>0, "Comparison requires nonempty candidate and model");
		stats.clear();stats.doTrace=true;final boolean candidateQuery=candidate.length<=model.length;
		final byte[] q=candidateQuery ? candidate : model, r=candidateQuery ? model : candidate;
		final float identity=QuantumAligner.alignAndTraceStatic(q, r, stats);
		require(Float.isFinite(identity) && identity>=0 && identity<=1 && Float.isFinite(stats.identity), "Quantum returned invalid identities");
		final int paired=stats.matches+stats.subs+stats.ns;
		require(paired+stats.ins==q.length && paired+stats.dels==stats.rStop-stats.rStart+1,
			"Trace operations must conserve whole query and reported reference span (QuantumAligner canonicalization contract)");
		require(stats.rStart>=0 && stats.rStop<r.length && paired<=candidate.length && paired<=model.length,
			"Reported pairing must fit both source sequences");
		return new Result(candidateQuery, identity, stats);
	}
	private static byte[] reference(String path, String id){
		final byte[][] result={null};stream(path, r->{if(first(r.id).equals(id)){require(result[0]==null, "Duplicate universal model: "+id);result[0]=bases(r.bases);}});
		require(result[0]!=null, "Missing exact universal model: "+id);return result[0];
	}
	private static byte[] bases(byte[] input){
		require(input!=null && input.length>0, "Empty sequence cannot be compared");
		final byte[] result=input.clone();
		for(int i=0; i<result.length; i++){byte b=result[i];if(b>='a' && b<='z'){b-=32;}require("ACGTURYSWKMBDHVN".indexOf(b)>=0, "Invalid sequence base");result[i]=b;}
		return result;
	}
	private static String first(String header){
		require(header!=null && !header.isEmpty(), "Input/reference identifier required");int end=0;
		while(end<header.length() && !Character.isWhitespace(header.charAt(end))){require(header.charAt(end)>32 && header.charAt(end)<127, "Locus key must be ASCII TSV-safe");end++;}
		require(end>0, "Empty first-token identifier");return header.substring(0, end);
	}
	private static void stream(String path, java.util.function.Consumer<Read> action){
		require(path!=null && new File(path).isFile(), "Missing FASTA: "+path);
		final Streamer st=StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 1);
		st.start();try{for(ListNum<Read> batch; (batch=st.nextList())!=null;){for(Read r:batch.list){action.accept(r);}st.returnList(batch);}require(!st.errorState(), "FASTA stream error: "+path);}
		finally{require(!ReadWrite.closeStream(st), "FASTA close error: "+path);}
	}
	private static ByteStreamWriter writer(String path){final ByteStreamWriter w=new ByteStreamWriter(FileFormat.testOutput(path, FileFormat.TEXT, null, true, false, false, false));w.start();return w;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	@Override public void accumulate(Worker w){if(w.failure!=null){w.failure.printStackTrace(outstream);}error|=!w.success;}
	@Override public boolean success(){return !error;}
	@Override public ReadWriteLock rwlock(){return lock;}
	static final class Worker extends Thread {
		Worker(NcrnaMitoSsuCompare owner_){owner=owner_;}
		@Override public void run(){try{
			final AlignmentStats stats=new AlignmentStats(true);
			for(int start; (start=owner.next.getAndAdd(16))<owner.unique.size();){
				final int end=Math.min(start+16, owner.unique.size());
				for(int i=start; i<end; i++){
					final Key key=owner.unique.get(i);for(int m=0; m<2; m++){key.results[m]=compare(key.bases, owner.refs[m], stats);}
				}
				final int done=owner.completed.addAndGet(end-start);
				if(done%1024==0 || done==owner.unique.size()){owner.outstream.println("Compared unique sequences: "+done+"/"+owner.unique.size());}
			}success=true;
		}catch(Throwable t){failure=t;}}
		final NcrnaMitoSsuCompare owner;boolean success;Throwable failure;
	}
	public static final class Result {
		Result(boolean cq, float id, AlignmentStats s){candidateQuery=cq;identity=id;traceIdentity=s.identity;start=s.rStart;stop=s.rStop;matches=s.matches;subs=s.subs;ins=s.ins;dels=s.dels;ns=s.ns;paired=matches+subs+ns;}
		public final boolean candidateQuery;public final float identity,traceIdentity;public final int start,stop,matches,subs,ins,dels,ns,paired;
	}
	private static final class Key {
		Key(byte[] b){bases=b;hash=Arrays.hashCode(b);}
		@Override public int hashCode(){return hash;}
		@Override public boolean equals(Object o){return o instanceof Key && Arrays.equals(bases, ((Key)o).bases);}
		final byte[] bases;final int hash;int copies;final Result[] results=new Result[2];
	}
	private static final class Member {Member(String n, int i){name=n;index=i;}final String name;final int index;}
	private String input,ref16,ref18,prefix;private int expected;private boolean error;private final PrintStream outstream;
	private final ObjectIntMap<Key> cache=new ObjectIntMap<Key>(Key.class);
	private final ArrayList<Key> unique=new ArrayList<Key>();private final ArrayList<Member> members=new ArrayList<Member>();
	private final byte[][] refs=new byte[2][];private final AtomicInteger next=new AtomicInteger(), completed=new AtomicInteger();private final ReadWriteLock lock=new ReentrantReadWriteLock();
	private static final String[] MODELS={"universal", "euk18S_universal"};
	private static final String[] SUFFIXES={".pairs.tsv.gz", ".members.tsv.gz", ".unique.fa.gz", ".summary.tsv"};
}
