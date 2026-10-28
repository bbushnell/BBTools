package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;

import consensus.BaseGraph;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import idaligner.ScrabbleAligner;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.FloatList;
import template.Accumulator;
import template.ThreadWaiter;

/** Checks consensus/HBM parity and exhaustive member-to-library identity.
 * Measures construction adequacy using the shorter query, not caller recall.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaResourceReport implements Accumulator<NcrnaResourceReport.ProcessThread> {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaResourceReport.class, false);
		String in=null, ref=null, hbm=null, assignments=null, out=null;int expected=-1, threads=1;float cut=Float.NaN;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("in")){in=value;}else if(key.equals("ref")){ref=value;}else if(key.equals("hbm")){hbm=value;}
			else if(key.equals("assignments")){assignments=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("expected")){expected=Parse.parseIntKMG(value);}else if(key.equals("reportcut")){cut=Float.parseFloat(value);}
			else if(key.equals("t") || key.equals("threads")){threads=Parse.parseIntKMG(value);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(in!=null && ref!=null && hbm!=null && assignments!=null && out!=null && expected>0 && cut>=0 && cut<=1,
			"Bind construction FASTA, consensus/HBM, builder membership, denominator and descriptive identity threshold");
		require(threads>0, "Reporter worker count must be positive");
		Shared.setThreads(threads);Shared.TRIM_READ_DESCRIPTION=false;
		new NcrnaResourceReport().run(in, ref, hbm, assignments, out, expected, cut);Shared.closeStream(pp.outstream);
	}

	private void run(final String input, final String refs, final String hbm, final String assignments, final String out, final int expected, final float cut){
		final ArrayList<Read> library=new ArrayList<Read>();RrnaResourceIO.stream(refs, library::add);
		final BaseGraph[] graphs=TrnaConsensusBuilder.loadModels(hbm);
		require(!library.isEmpty() && graphs!=null && library.size()==graphs.length, "Consensus and HBM counts must agree before reporting coverage");
		final ObjectIntMap<String> modelIds=new ObjectIntMap<String>(String.class);
		for(int i=0; i<graphs.length; i++){
			final Read r=library.get(i);final String id=RrnaResourceIO.first(r.id);
			require(!modelIds.contains(id) && r.bases!=null && r.length()>0 && id.equals(RrnaResourceIO.first(graphs[i].name))
				&& Arrays.equals(r.bases, graphs[i].original), "Model names/order/reference bases must match exactly: "+id);
			modelIds.put(id, i);
		}
		final ObjectMap<String,Member> members=load(assignments, modelIds, expected);
		RrnaResourceIO.fresh(out, ".members.tsv.gz", ".summary.tsv", ".models.tsv", ".histogram.tsv");
		final ArrayList<Read> queries=new ArrayList<Read>(expected);
		final Member[] ordered=new Member[expected];
		RrnaResourceIO.stream(input, r->{
			final String id=RrnaResourceIO.first(r.id);final Member m=members.get(id);
			require(m!=null && !m.seen && r.id.equals(m.header) && r.bases!=null && r.length()>0,
				"Each original builder member must occur once with unchanged header and nonempty sequence: "+id);
			require(queries.size()<expected, "Input contains more records than the declared assignment denominator");
			m.seen=true;ordered[queries.size()]=m;queries.add(r);
		});
		require(queries.size()==expected, "Input is missing original builder members: "+queries.size()+" expected="+expected);
		final float[] bestIdentity=new float[expected], ownIdentity=new float[expected];final int[] winners=new int[expected];
		final AtomicInteger next=new AtomicInteger();
		final int threads=Math.min(Shared.threads(), 1+(expected-1)/BATCH_SIZE);
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){workers.add(new ProcessThread(queries, library, ordered, next, bestIdentity, ownIdentity, winners));}
		final boolean completed=ThreadWaiter.startAndWait(workers, this);
		if(!completed){throw new IllegalStateException("Member alignment failed; no report was written", failure);}
		require(alignmentCalls==(long)expected*library.size(), "Every member must be aligned to every model before writing coverage");
		final ByteStreamWriter rows=RrnaResourceIO.writer(out+".members.tsv.gz");final ByteBuilder b=new ByteBuilder();
		final FloatList values=new FloatList();final long[] histogram=new long[1001];
		final int[] assigned=new int[library.size()], bestCounts=new int[library.size()], states=new int[4];
		try{
			rows.println("locus_id\tlength\tbuilder_status\tassigned_model\tassigned_identity\tbest_model\tbest_shorter_query_identity\tat_report_cut");
			for(int q=0; q<expected; q++){
				final Read r=queries.get(q);final String id=RrnaResourceIO.first(r.id);final Member m=ordered[q];
				final float best=bestIdentity[q], own=ownIdentity[q];final int winner=winners[q];
				if(m.model>=0){assigned[m.model]++;states[0]++;}else{states[m.status.equals("small")?1:2]++;}
				if(best>=cut){states[3]++;}bestCounts[winner]++;values.add(best);histogram[Math.min(1000, (int)(1000*best))]++;
				b.clear().append(id).tab().append(r.length()).tab().append(m.status).tab();
				if(m.model>=0){b.append(RrnaResourceIO.first(library.get(m.model).id)).tab().append(own, 6);}
				else{b.append("NA\tNA");}
				rows.print(b.tab().append(RrnaResourceIO.first(library.get(winner).id)).tab().append(best, 6).tab().append(best>=cut).nl());
			}
		}finally{RrnaResourceIO.finish(rows);}
		require(values.size==expected && states[0]+states[1]+states[2]==expected && alignmentCalls==(long)expected*library.size(),
			"Exhaustive identity and builder disposition must conserve every input member, including orphans");
		values.sort();final ByteStreamWriter summary=RrnaResourceIO.writer(out+".summary.tsv"), models=RrnaResourceIO.writer(out+".models.tsv"), hist=RrnaResourceIO.writer(out+".histogram.tsv");
		try{
			summary.println("members\tmodels\tkept\tsmall\torphan\treport_cut\tcovered\tuncovered\tmin\tp10\tmedian\tp90\tmax\talignment_calls\tinterpretation");
			b.clear().append(expected).tab().append(library.size()).tab().append(states[0]).tab().append(states[1]).tab().append(states[2])
				.tab().append(cut, 4).tab().append(states[3]).tab().append(expected-states[3]);
			for(int percentile:new int[]{0, 10, 50, 90, 100}){b.tab().append(values.array[(int)((expected-1L)*percentile/100)], 6);}
			summary.print(b.tab().append(alignmentCalls).append("\tIN_CORPUS_SHORTER_QUERY_IDENTITY_NOT_CALLER_RECALL\n"));
			models.println("model_id\tlength\tbuilder_assigned\tbest_match_members");
			for(int i=0; i<library.size(); i++){models.print(b.clear().append(RrnaResourceIO.first(library.get(i).id)).tab().append(library.get(i).length()).tab().append(assigned[i]).tab().append(bestCounts[i]).nl());}
			hist.println("identity_bin_lower\tmembers");
			for(int i=0; i<histogram.length; i++){if(histogram[i]>0){hist.print(b.clear().append(i/1000.0, 3).tab().append(histogram[i]).nl());}}
		}finally{RrnaResourceIO.finish(summary, models, hist);}
		System.err.println("NCRNA_RESOURCE_REPORT_PASS members="+expected+" models="+library.size()+" covered="+states[3]+"; exhaustive in-corpus identity, not caller recall");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void accumulate(final ProcessThread worker){
		assert(!worker.isAlive()):"ThreadWaiter must join before publishing member identities";
		alignmentCalls+=worker.calls;
		errorState|=!worker.success;
		if(failure==null){failure=worker.failure;}
	}
	@Override
	public boolean success(){return !errorState;}
	@Override
	public ReadWriteLock rwlock(){return null;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Each worker owns complete query rows; publication waits for every worker. */
	static final class ProcessThread extends Thread {
		ProcessThread(final ArrayList<Read> q, final ArrayList<Read> r, final Member[] m,
			final AtomicInteger n, final float[] b, final float[] o, final int[] w){
			queries=q;library=r;members=m;next=n;bestIdentity=b;ownIdentity=o;winners=w;
			assert(q.size()==m.length && b.length==m.length && o.length==m.length && w.length==m.length):
				"Disjoint result slots must cover exactly the validated input members";
		}
		@Override
		public void run(){
			try{
				for(int first; (first=next.getAndAdd(BATCH_SIZE))<queries.size(); ){
					final int end=Math.min(queries.size(), first+BATCH_SIZE);
					for(int q=first; q<end; q++){align(q);}
				}
				success=true;
			}catch(Throwable t){
				failure=t;next.set(queries.size());
			}
		}
		private void align(final int q){
			final Read read=queries.get(q);final Member member=members[q];
			float best=-1, own=-1;int winner=-1;
			for(int i=0; i<library.size(); i++){
				final byte[] ref=library.get(i).bases;
				final float identity=ScrabbleAligner.alignStatic(read.bases, ref, null);
				if(!Float.isFinite(identity) || identity<0 || identity>1){throw new IllegalArgumentException("Nonfractional alignment identity: "+read.id);}
				calls++;
				if(identity>best){best=identity;winner=i;}if(i==member.model){own=identity;}
			}
			assert(winner>=0):"A nonempty library must yield a diagnostic model for every member";
			bestIdentity[q]=best;ownIdentity[q]=own;winners[q]=winner;
		}
		private final ArrayList<Read> queries, library;
		private final Member[] members;
		private final AtomicInteger next;
		private final float[] bestIdentity, ownIdentity;
		private final int[] winners;
		private long calls;
		private boolean success;
		private Throwable failure;
	}

	private static ObjectMap<String,Member> load(final String path, final ObjectIntMap<String> models, final int expected){
		final ObjectMap<String,Member> members=new ObjectMap<String,Member>(String.class, Member.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==6, "Builder assignment schema changed");
				final String header=p.parseString(0), id=RrnaResourceIO.first(header), status=p.parseString(5), label=p.parseString(4);
				require(!members.contains(id) && (status.equals("kept") || status.equals("small") || status.equals("orphan")), "Members need unique identities and explicit dispositions: "+id);
				int model=-1;
				if(status.equals("kept")){
					final String key=RrnaResourceIO.first(label);require(models.contains(key), "Kept assignment names a missing model: "+key);model=models.get(key);
				}else{require(label.isEmpty(), "Non-kept members must not claim a model assignment: "+id);}
				members.put(id, new Member(header, status, model));
			}
		}finally{require(!in.close(), "Builder assignment read failed");}
		require(members.size()==expected, "Assignment denominator changed: "+members.size()+" expected="+expected);return members;
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Member {
		Member(final String h, final String s, final int m){header=h;status=s;model=m;}
		final String header, status;final int model;boolean seen;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private long alignmentCalls;
	private boolean errorState;
	private Throwable failure;
	private static final int BATCH_SIZE=32;
}
