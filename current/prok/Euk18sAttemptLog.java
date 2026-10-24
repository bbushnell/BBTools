package prok;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import idaligner.AlignmentStats;
import structures.ByteBuilder;

/** Diagnostic attempts from the production 18S alignment, without secondary alignment.
 * @author Ganyu */
public final class Euk18sAttemptLog {
	static final String PROPERTY="bbtools.diagnostic.euk18sAttemptsOut";
	private static final String OUTPUT=System.getProperty(PROPERTY);
	static final String HEADER="family\tcontig\tstrand\tpass\twindowStart\twindowStop\tcontigLength\troute\tclipModel\tattemptOrdinal\tmodel\tmodelName\tshared9\tbestShared9\tcountRank\tidentity\talignedStart\talignedStop\talignedLength\tvalidBounds\tidCutoff\tlengthAllowed\taccepted";
	public static Euk18sAttemptLog forFamily(String family){
		return OUTPUT!=null && Euk18sRuntimeConfig.FAMILY.equals(family) ? Holder.INSTANCE : null;
	}
	private static final class Holder{
		static final Euk18sAttemptLog INSTANCE=open();
		static Euk18sAttemptLog open(){
			final Euk18sAttemptLog result=new Euk18sAttemptLog(OUTPUT);
			Runtime.getRuntime().addShutdownHook(new Thread("18S-attempt-log-exit"){
				@Override public void run(){try{result.close();}catch(Throwable e){e.printStackTrace(System.err);Runtime.getRuntime().halt(1);}}
			});
			return result;
		}
	}
	public Euk18sAttemptLog(String path){
		require(path!=null && !path.isEmpty(), "Attempt log needs a fresh plain TSV path");
		RrnaResourceIO.fresh(path, "");
		try{
			// Synchronous buffered output: no live writer thread can prevent the shutdown hook.
			out=new BufferedOutputStream(Files.newOutputStream(Paths.get(path),StandardOpenOption.CREATE_NEW,StandardOpenOption.WRITE));
			write("# Coordinates are zero-based inclusive on the strand-oriented contig; raw invalid alignment bounds are retained with validBounds=false. Missing alignment is NA.\n");
			write("# shared9 uses unique query keys weighted by model occurrence multiplicity; countRank is across all models, count descending/model ID ascending, including zeros.\n");
			write("# accepted means per-attempt identity and length eligibility, before model selection, endpoint refinement, overlap resolution or final path selection; it is not final GFF retention.\n");
			write("# attemptOrdinal resets for each window/clipModel verifier phase; joined forces a model and does not imply countRank=1. Diagnostic-only joined score refresh adds mapping work, never alignment.\n");
			write(HEADER+"\n");out.flush();
		}catch(IOException e){throw new IllegalStateException("Cannot create 18S attempt log: "+path,e);}
	}
	synchronized void attempt(String family,String contig,int strand,int pass,int ws,int we,int contigLength,
			boolean joined,boolean clip,int ordinal,int model,String modelName,TrnaKmerIndex index,int modelCount,
			float identity,AlignmentStats stats,boolean validBounds,float cutoff,boolean lengthAllowed){
		require(!closed && Euk18sRuntimeConfig.FAMILY.equals(family) && strand>=0 && strand<=1 && pass>=1
			&& ws>=0 && we>=ws && we<contigLength && ordinal>0 && model>=0 && model<modelCount && index!=null,
			"Attempt must retain an actual 18S model alignment and physical candidate window");
		require(Float.isFinite(identity) && Float.isFinite(cutoff) && (!validBounds || stats!=null),
			"Attempt requires actual finite identity/cutoff; valid coordinates require a trace");
		cell(contig);cell(modelName);
		final int count=index.lastSharedCount(model), best=index.lastMaxShared(), rank=countRank(index,model,modelCount);
		require(count>=0 && best>=count, "Shared scores must belong to the current candidate window");
		b.clear().append(family).tab().append(contig).tab().append(strand==0?'+':'-').tab().append(pass)
			.tab().append(ws).tab().append(we).tab().append(contigLength).tab().append(joined?"joined":"ordinary")
			.tab().append(clip).tab().append(ordinal).tab().append(model).tab().append(modelName)
			.tab().append(count).tab().append(best).tab().append(rank).tab().append(Float.toString(identity)).tab();
		if(stats==null){b.append("NA\tNA\tNA");}
		else{b.append(ws+stats.rStart).tab().append(ws+stats.rStop).tab().append(stats.rStop-stats.rStart+1);}
		b.tab().append(validBounds).tab().append(Float.toString(cutoff)).tab().append(lengthAllowed)
			.tab().append(identity>=cutoff && lengthAllowed).nl();
		try{out.write(b.array,0,b.length());rows++;}
		catch(IOException e){throw new IllegalStateException("18S attempt row write failed",e);}
	}
	public static int countRank(TrnaKmerIndex index,int model,int models){
		require(index!=null && models>0 && model>=0 && model<models,"Rank needs a member of the current model library");
		final int count=index.lastSharedCount(model);int rank=1;
		for(int i=0;i<models;i++){final int other=index.lastSharedCount(i);if(other>count || (other==count && i<model)){rank++;}}
		return rank;
	}
	public synchronized void close(){
		require(!closed,"Attempt writer must close once after all calling threads finish");
		try{out.close();closed=true;}catch(IOException e){throw new IllegalStateException("18S attempt log close failed",e);}
	}
	public synchronized long rows(){return rows;}
	private void write(String text)throws IOException{out.write(text.getBytes(StandardCharsets.UTF_8));}
	private static void cell(String value){require(value!=null && value.indexOf('\t')<0 && value.indexOf('\n')<0 && value.indexOf('\r')<0,"Attempt labels must fit one TSV cell");}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final OutputStream out;
	private final ByteBuilder b=new ByteBuilder(512);
	private boolean closed;
	private long rows;
}
