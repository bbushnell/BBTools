package prok;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** One compact row per window reaching the mitochondrial tier after original
 * verification fails. Observation never creates an alignment or changes ranking.
 * @author Brian, Raiden */
public final class MitoSsuTierLog {
	static synchronized MitoSsuTierLog configured(){
		final String path=System.getProperty("bbtools.diagnostic.mitoTierOut");
		if(path==null){return null;}
		if(instance==null){instance=new MitoSsuTierLog(path);}
		return instance;
	}
	public MitoSsuTierLog(final String path){
		if(path==null || path.isEmpty()){throw new IllegalArgumentException("Tier observation requires a fresh TSV path");}
		RrnaResourceIO.fresh(path,"");out=RrnaResourceIO.writer(path);out.println(HEADER);
	}
	public synchronized void record(final String contig, final int strand, final int pass, final int start, final int stop,
			final int queryUnique, final int top, final int second, final int third, final int models,
			final long alignments, final long skips, final int acceptedModel, final int acceptedRank,
			final int acceptedAttempt, final int acceptedShared, final String outcome, final long firstTwoCalls){
		if(closed || contig==null || contig.indexOf('\t')>=0 || contig.indexOf('\n')>=0 || contig.indexOf('\r')>=0
			|| strand<0 || strand>1 || pass<1 || start<0 || stop<start || queryUnique<0 || top<second || second<third || third<0
			|| models<0 || alignments<models || skips<0 || (acceptedModel>=0)!=(acceptedAttempt>0)
			|| firstTwoCalls<0 || firstTwoCalls>alignments || acceptedAttempt>models || (acceptedModel>=0 && (acceptedRank<1 || acceptedShared<1))){
			throw new IllegalArgumentException("Tier row must bind real window scores, actual model attempts and accepted-model provenance");
		}
		b.clear().append(contig).tab().append(strand).tab().append(pass).tab().append(start).tab().append(stop)
			.tab().append(queryUnique).tab().append(top).tab().append(second).tab().append(third).tab().append(models)
			.tab().append(alignments).tab().append(skips).tab().append(acceptedModel).tab().append(acceptedRank)
			.tab().append(acceptedAttempt).tab().append(acceptedShared).tab().append(outcome).tab().append(firstTwoCalls).nl();out.print(b);
	}
	public synchronized void close(){
		if(closed){throw new IllegalStateException("Tier observer must close once after all caller workers join");}
		RrnaResourceIO.finish(out);closed=true;
	}
	public static synchronized void closeConfigured(){if(instance!=null){instance.close();instance=null;}}
	public static final String HEADER="contig\tstrand\tpass\twindow_start\twindow_stop\tquery_unique9\ttop_shared9\tsecond_shared9\tthird_shared9\tmodels_attempted\taligner_calls\tlength_skips\taccepted_model\taccepted_rank\taccepted_attempt\taccepted_shared9\toutcome\tfirst_two_aligner_calls";
	private static MitoSsuTierLog instance;
	private final ByteStreamWriter out;
	private final ByteBuilder b=new ByteBuilder(256);
	private boolean closed;
}
