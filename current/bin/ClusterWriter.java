package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

public class ClusterWriter{
	
	ClusterWriter(PrintStream outstream_, boolean writeChaff_, boolean loud_, boolean overwrite_, boolean append_){
		outstream=outstream_;
		overwrite=overwrite_;
		append=append_;
		writeChaff=writeChaff_;
		loud=loud_;
	}

	/**
	 * Writes clustered contigs to output files using specified pattern.
	 * Supports individual bin files with % substitution or single combined output.
	 * Optionally writes small clusters to a separate chaff file.
	 *
	 * @param pattern Output file pattern with % for bin numbering
	 * @param clusters List of bins/clusters to output
	 * @param minBases Minimum base count threshold for main output
	 * @param minContigs Minimum contig count threshold for main output
	 */
	void outputClusters(String pattern, String sizeHist, ArrayList<? extends Bin> clusters,
			long minBases, int minContigs, Binner binner, DataLoader loader) {
		Timer t=new Timer();
		
		if(pattern!=null) {
			String name=ReadWrite.stripPath(pattern);
			final boolean isFile=FileFormat.hasFastaOrFastqExtension(name) || 
				name.contains("%") || name.contains("#");
			if(isFile && !name.contains("%") && name.contains("#")) {
				pattern=Tools.replaceLastInstanceOf(pattern, '#', '%');//A common mistake I make
				//Maybe prudent to make sure it is in a filename rather than dirname
			}
			if(!isFile) {
				if(!pattern.endsWith("/")) {pattern=pattern+"/";}
				pattern=pattern+"bin_%.fa";
			}
			outstream.println("Writing clusters to "+pattern);
		}
		long sizeOverLimit=0;
		if(pattern!=null && pattern.indexOf('%')>=0) {
			
			ByteStreamWriter chaff=null;
			if(writeChaff) {
				chaff=ByteStreamWriter.makeBSW(pattern.replaceFirst("%", "chaff"), overwrite, append, true);
			}
			
			final ByteBuilder bb=new ByteBuilder(32768);
			for(int i=0; i<clusters.size(); i++) {
				Bin a=clusters.get(i);
				if(a.size()>=minBases && a.numContigs()>=minContigs) {
					String fname=pattern.replaceFirst("%", Integer.toString(i));
					if(fname.contains("%contam")) {fname=fname.replaceFirst("%contam", String.format("con%.4f", a.contam));}
					else if(fname.contains("contam%")) {fname=fname.replaceFirst("contam%", String.format("con%.4f", a.contam));}
					if(fname.contains("%comp")) {fname=fname.replaceFirst("%comp", String.format("com%.4f", a.completeness));}
					else if(fname.contains("comp%")) {fname=fname.replaceFirst("comp%", String.format("com%.4f", a.completeness));}
					final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(fname, overwrite, append, a.size()>=MIN_SIZE_FOR_SUBPROCESS);
					printBin(a, bsw, bb, -1);
					bsw.poison();
					clustersWritten++;
					contigsWritten+=a.numContigs();
					basesWritten+=a.size();
				}else {
					printBin(a, chaff, bb, i+1);
				}
			}
			if(chaff!=null) {chaff.poisonAndWait();}
		}else {
			final ByteBuilder bb=new ByteBuilder(32768);
			final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(pattern, overwrite, append, true);
			for(int i=0; i<clusters.size(); i++) {
				Bin a=clusters.get(i);
				printBin(a, bsw, bb, i+1);
				clustersWritten++;
				contigsWritten+=a.numContigs();
				basesWritten+=a.size();
			}
			if(bsw!=null) {bsw.poisonAndWait();}
		}
		float cpct=contigsWritten*100f/loader.contigsLoaded;
		float bpct=basesWritten*100f/loader.basesLoaded;
		if(loud) {
			outstream.println("\nMetric   \t        In\t       Out\tPercent");
			outstream.println("Clusters \t"+Tools.padLeft(0, 10)+"\t"+Tools.padLeft(clustersWritten, 10));
			outstream.println("Contigs  \t"+Tools.padLeft(loader.contigsLoaded, 10)+"\t"+
					Tools.padLeft(contigsWritten, 10)+"\t"+Tools.format("%.2f%%", cpct));
			outstream.println("Bases    \t"+Tools.padLeft(loader.basesLoaded, 10)+"\t"+
					Tools.padLeft(basesWritten, 10)+"\t"+Tools.format("%.2f%%", bpct));
		}
		if(sizeHist!=null) {
			final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(sizeHist, overwrite, append, true);
			bsw.print("#Size\tGood\tBad\n");
			for(int i=0; i<binner.goodMergeSize.length; i++) {
				long size=1L<<i;
				bsw.print(size).tab().print(binner.goodMergeSize[i]).tab().print(binner.badMergeSize[i]).nl();
			}
			bsw.poisonAndWait();
		}
		t.stop("Time:");
	}
	
	/**
	 * Prints a bin to output stream, handling both clusters and single contigs.
	 * Dispatches to appropriate print method based on bin type.
	 *
	 * @param a Bin to print (cluster or single contig)
	 * @param bsw Output stream writer (may be null)
	 * @param bb Byte buffer for efficient output
	 * @param id Numeric identifier for the bin
	 */
	private void printBin(Bin a, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		if(a.isCluster()) {printCluster((Cluster)a, bsw, bb, id);}
		else {printContig((Contig)a, bsw, bb, id);}
	}
	
	/**
	 * Prints all contigs in a cluster to the output stream.
	 * Sorts contigs within cluster and buffers output for efficiency.
	 *
	 * @param a Cluster containing multiple contigs
	 * @param bsw Output stream writer (may be null)
	 * @param bb Byte buffer for batched output
	 * @param id Numeric identifier for the cluster
	 */
	private void printCluster(Cluster a, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		ArrayList<Contig> contigs=a.contigs;
		Collections.sort(contigs);
		for(Contig c : contigs) {
			c.appendTo(bb, id);
			if(bb.length>=65536) {
				if(bsw!=null) {bsw.print(bb);}
				bb.clear();
			}
		}
		if(bsw!=null && !bb.isEmpty()) {bsw.print(bb);}
		bb.clear();
	}
	
	/**
	 * Prints a single contig to the output stream.
	 *
	 * @param c Contig to print
	 * @param bsw Output stream writer (may be null)
	 * @param bb Byte buffer for output formatting
	 * @param id Numeric identifier for the contig
	 */
	private void printContig(Contig c, ByteStreamWriter bsw, ByteBuilder bb, int id) {
		c.appendTo(bb, id);
		if(bsw!=null && !bb.isEmpty()) {bsw.print(bb);}
		bb.clear();
	}
	
	private PrintStream outstream;
	
	long clustersWritten=0;
	long contigsWritten=0;
	long basesWritten=0;
	final boolean overwrite;
	final boolean append;
	final boolean writeChaff;
	final boolean loud;
	
	static int MIN_SIZE_FOR_SUBPROCESS=1000000;
	
}
