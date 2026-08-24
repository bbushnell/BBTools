package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.GlocalAligner;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.Read;
import structures.ListNum;

/**
 * Calibration and coverage measurement tool for ncRNA consensus libraries.
 * A mini CallGenes for calibrating each individual ncRNA class.
 *
 * For each query, shortlists reference models via a 7-mer inverted index
 * (TrnaKmerIndex), aligns in decreasing kmer-sharing order, and reports
 * mapping/alignment pass rates plus diagnostic histograms.
 *
 * Output follows BBDuk conventions: out/outu= captures queries that FAIL
 * (unmapped or below identity threshold), outm= captures queries that PASS.
 *
 * @author Noire, Brian Bushnell
 */
public class NcrnaCalibrator {

	public static void main(String[] args){
		Timer t=new Timer();
		NcrnaCalibrator x=new NcrnaCalibrator(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public NcrnaCalibrator(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("ref") || a.equals("library")){
				refFile=b;
			}else if(a.equals("outm") || a.equals("outcovered")){
				outMatched=b;
			}else if(a.equals("trim")){
				trim=Integer.parseInt(b);
			}else if(a.equals("minid") || a.equals("minidentity")){
				minId=Float.parseFloat(b);
			}else if(a.equals("minkmers") || a.equals("minkmerhits")){
				minKmers=Integer.parseInt(b);
			}else if(a.equals("aligner")){
				alignerType=b.toLowerCase();
			}else if(a.equals("indexk") || a.equals("ik")){
				indexK=Integer.parseInt(b);
			}else if(a.equals("topn") || a.equals("indextopn")){
				indexTopN=Integer.parseInt(b);
			}else if(a.equals("adaptive")){
				adaptive=Parse.parseBoolean(b);
			}else if(a.equals("floor") || a.equals("adaptfloor")){
				adaptFloor=Float.parseFloat(b);
			}else if(a.equals("topfrac") || a.equals("adapttopfrac")){
				adaptTopFrac=Float.parseFloat(b);
			}else if(a.equals("qfrac") || a.equals("adaptqfrac")){
				adaptQFrac=Float.parseFloat(b);
			}else if(a.equals("fixedminhits")){
				fixedMinHits=Integer.parseInt(b);
			}else if(a.equals("anihist")){
				aniHistFile=b;
			}else if(a.equals("kmerhist")){
				kmerHistFile=b;
			}else if(a.equals("ratiohist")){
				ratioHistFile=b;
			}else if(a.equals("histbins")){
				histBins=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//handled
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		in=parser.in1;
		out=parser.out1;
		if(in==null){throw new RuntimeException("Error - in= is required.");}
		if(refFile==null){throw new RuntimeException("Error - ref= is required.");}
		overwrite=parser.overwrite;

		ffin=FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
		ffout=(out!=null ? FileFormat.testOutput(out, FileFormat.FASTA, null, true, overwrite, false, false) : null);
		ffoutm=(outMatched!=null ? FileFormat.testOutput(outMatched, FileFormat.FASTA, null, true, overwrite, false, false) : null);
	}

	void process(Timer t){
		byte[][] library=loadLibrary(refFile);
		final int nModels=library.length;
		outstream.println("Loaded "+nModels+" reference models from "+refFile);

		TrnaKmerIndex kmerIndex=new TrnaKmerIndex(library, indexK, adaptive,
				adaptFloor, adaptTopFrac, adaptQFrac, fixedMinHits);

		ConcurrentReadOutputStream rosFailed=null, rosPassed=null;
		if(ffout!=null){
			rosFailed=ConcurrentReadOutputStream.getStream(ffout, null, 4, null, false);
			rosFailed.start();
		}
		if(ffoutm!=null){
			rosPassed=ConcurrentReadOutputStream.getStream(ffoutm, null, 4, null, false);
			rosPassed.start();
		}

		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ffin, null);
		cris.start();

		long totalQueries=0, passedMapping=0, failedMapping=0, passedAlignment=0, failedAlignment=0;
		long totalAlignments=0;

		final int aniBins=histBins;
		final long[] aniHist=new long[aniBins+1];
		final long[] topKmerHist=new long[1001];
		final long[] acceptedKmerHist=new long[1001];
		final long[] ratioHist=new long[101];

		ArrayList<Read> failedBatch=new ArrayList<>();
		ArrayList<Read> passedBatch=new ArrayList<>();
		long listId=0;

		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases==null || r.length()==0){continue;}
				totalQueries++;

				byte[] seq=r.bases;
				byte[] outSeq=(trim>0 && seq.length>2*trim ?
						Arrays.copyOfRange(seq, trim, seq.length-trim) : seq);

				int[] shortlist=kmerIndex.shortlist(seq, indexTopN);
				int topKmerCount=(shortlist.length>0 ? kmerIndex.lastSharedCount(shortlist[0]) : 0);

				if(shortlist.length==0 || topKmerCount<minKmers){
					failedMapping++;
					addToAniHist(aniHist, 0, aniBins);
					if(rosFailed!=null){
						failedBatch.add(new Read(outSeq, null, r.id, r.numericID));
					}
					continue;
				}
				passedMapping++;
				topKmerHist[Tools.min(topKmerCount, topKmerHist.length-1)]++;

				int[] allKmerCounts=new int[shortlist.length];
				for(int j=0; j<shortlist.length; j++){
					allKmerCounts[j]=kmerIndex.lastSharedCount(shortlist[j]);
				}

				float bestId=0;
				int acceptedModel=-1;
				int acceptedKmerCount=0;

				for(int j=0; j<shortlist.length; j++){
					int m=shortlist[j];
					float id=align(library[m], seq);
					totalAlignments++;
					if(id>bestId){bestId=id;}
					if(id>=minId && acceptedModel<0){
						acceptedModel=m;
						acceptedKmerCount=allKmerCounts[j];
						break;
					}
				}

				addToAniHist(aniHist, bestId, aniBins);

				if(acceptedModel>=0){
					passedAlignment++;
					acceptedKmerHist[Tools.min(acceptedKmerCount, acceptedKmerHist.length-1)]++;
					int ratio100=(topKmerCount>0 ? (int)Math.round(100.0*acceptedKmerCount/topKmerCount) : 100);
					ratioHist[Tools.min(ratio100, 100)]++;
					if(rosPassed!=null){
						passedBatch.add(new Read(outSeq, null, r.id, r.numericID));
					}
				}else{
					failedAlignment++;
					if(rosFailed!=null){
						failedBatch.add(new Read(outSeq, null, r.id, r.numericID));
					}
				}

				if(verbose && totalQueries%10000==0){
					outstream.println("  Processed "+totalQueries+", passed="+passedAlignment+
							" ("+String.format("%.2f", 100.0*passedAlignment/totalQueries)+"%)");
				}
			}
			cris.returnList(ln);

			if(failedBatch.size()>=1000){
				if(rosFailed!=null){rosFailed.add(new ArrayList<>(failedBatch), listId++);}
				failedBatch.clear();
			}
			if(passedBatch.size()>=1000){
				if(rosPassed!=null){rosPassed.add(new ArrayList<>(passedBatch), listId++);}
				passedBatch.clear();
			}
		}
		ReadWrite.closeStream(cris);

		if(!failedBatch.isEmpty() && rosFailed!=null){rosFailed.add(failedBatch, listId++);}
		if(!passedBatch.isEmpty() && rosPassed!=null){rosPassed.add(passedBatch, listId++);}
		if(rosFailed!=null){errorState|=ReadWrite.closeStream(rosFailed);}
		if(rosPassed!=null){errorState|=ReadWrite.closeStream(rosPassed);}

		outstream.println();
		outstream.println("Passed mapping:   \t"+passedMapping+"\t"+
				String.format("%.3f", 100.0*passedMapping/Math.max(1, totalQueries))+"%");
		outstream.println("Failed mapping:   \t"+failedMapping+"\t"+
				String.format("%.3f", 100.0*failedMapping/Math.max(1, totalQueries))+"%");
		outstream.println("Passed alignment: \t"+passedAlignment+"\t"+
				String.format("%.3f", 100.0*passedAlignment/Math.max(1, passedMapping))+"%");
		outstream.println("Failed alignment: \t"+failedAlignment+"\t"+
				String.format("%.3f", 100.0*failedAlignment/Math.max(1, passedMapping))+"%");
		outstream.println();
		outstream.println("Queries:          \t"+totalQueries);
		outstream.println("Total alignments: \t"+totalAlignments);
		outstream.println("Models:           \t"+nModels);
		outstream.println("Aligner:          \t"+alignerType);
		outstream.println("Index k:          \t"+indexK);
		outstream.println("Min kmers:        \t"+minKmers);
		outstream.println("Min identity:     \t"+String.format("%.4f", minId));
		outstream.println("Trim:             \t"+trim);
		if(adaptive){
			outstream.println("Adaptive:         \tfloor="+adaptFloor+" topFrac="+adaptTopFrac+" qFrac="+adaptQFrac);
		}else{
			outstream.println("Fixed minHits:    \t"+fixedMinHits);
		}

		writeHist(buildAniHistText(aniHist, aniBins, totalQueries), aniHistFile, "ANI Histogram");
		writeHist(buildKmerHistText(topKmerHist, "TopKmers", passedMapping), kmerHistFile, "Top-Hit Shared Kmers Histogram");
		if(passedAlignment>0){
			writeHist(buildKmerHistText(acceptedKmerHist, "AcceptedKmers", passedAlignment),
					kmerHistFile!=null ? kmerHistFile.replace(".txt", "_accepted.txt").replace(".tsv", "_accepted.tsv") : null,
					"Accepted-Hit Shared Kmers Histogram");
			writeHist(buildRatioHistText(ratioHist, passedAlignment), ratioHistFile, "Kmer Ratio Histogram");
		}

		t.stop();
		outstream.println();
		outstream.println("Time:             \t"+t);
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state.");
		}
	}

	private float align(byte[] model, byte[] query){
		if(alignerType.equals("quantum") || alignerType.equals("q")){
			return QuantumAligner.alignStatic(model, query, null);
		}else if(alignerType.equals("scrabble") || alignerType.equals("s")){
			return ScrabbleAligner.alignStatic(model, query, null);
		}else if(alignerType.equals("glocal") || alignerType.equals("g")){
			return GlocalAligner.alignStatic(model, query, null);
		}else{
			throw new RuntimeException("Unknown aligner: "+alignerType+
					" (use quantum, scrabble, or glocal)");
		}
	}

	private byte[][] loadLibrary(String fname){
		ArrayList<byte[]> list=new ArrayList<>();
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases!=null && r.length()>0){
					Tools.toUpperCase(r.bases);
					list.add(r.bases);
				}
			}
			cris.returnList(ln);
		}
		ReadWrite.closeStream(cris);
		return list.toArray(new byte[0][]);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Histogram Output      ----------------*/
	/*--------------------------------------------------------------*/

	private void writeHist(String text, String file, String label){
		if(file!=null && file.equalsIgnoreCase("stderr")){
			outstream.println();
			outstream.println("=== "+label+" ===");
			outstream.print(text);
		}else if(file!=null && !file.equalsIgnoreCase("null")){
			ReadWrite.writeString(text, file);
		}else{
			outstream.println();
			outstream.println("=== "+label+" ===");
			outstream.print(text);
		}
	}

	private static void addToAniHist(long[] hist, float id, int bins){
		int bin=(int)(id*bins);
		hist[Tools.min(bin, bins)]++;
	}

	private String buildAniHistText(long[] hist, int bins, long total){
		StringBuilder sb=new StringBuilder();
		sb.append("#ANI\tCount\tCumulative%\n");
		long cum=0;
		for(int i=0; i<=bins; i++){
			if(hist[i]>0){
				cum+=hist[i];
				sb.append(String.format("%.3f", (float)i/bins)).append('\t')
					.append(hist[i]).append('\t')
					.append(String.format("%.3f", 100.0*cum/Math.max(1, total))).append("%\n");
			}
		}
		return sb.toString();
	}

	private String buildKmerHistText(long[] hist, String label, long total){
		StringBuilder sb=new StringBuilder();
		sb.append("#").append(label).append("\tCount\tCumulative%\n");
		long cum=0;
		int lastNonzero=0;
		for(int i=0; i<hist.length; i++){if(hist[i]>0){lastNonzero=i;}}
		for(int i=0; i<=lastNonzero; i++){
			cum+=hist[i];
			if(hist[i]>0){
				sb.append(i).append('\t').append(hist[i]).append('\t')
					.append(String.format("%.3f", 100.0*cum/Math.max(1, total))).append("%\n");
			}
		}
		return sb.toString();
	}

	private String buildRatioHistText(long[] hist, long total){
		StringBuilder sb=new StringBuilder();
		sb.append("#Ratio%\tCount\tCumulative%\n");
		long cum=0;
		for(int i=0; i<=100; i++){
			cum+=hist[i];
			if(hist[i]>0){
				sb.append(i).append('\t').append(hist[i]).append('\t')
					.append(String.format("%.3f", 100.0*cum/Math.max(1, total))).append("%\n");
			}
		}
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in;
	private String refFile;
	private String out;
	private String outMatched;
	private String aniHistFile;
	private String kmerHistFile;
	private String ratioHistFile;
	private final FileFormat ffin;
	private final FileFormat ffout;
	private final FileFormat ffoutm;
	private boolean overwrite=true;
	private boolean errorState=false;

	private int trim=0;
	private float minId=0.75f;
	private int minKmers=0;
	private String alignerType="quantum";
	private int indexK=7;
	private int indexTopN=60;
	private boolean adaptive=true;
	private float adaptFloor=11;
	private float adaptTopFrac=0.48f;
	private float adaptQFrac=0.072f;
	private int fixedMinHits=12;
	private int histBins=100;

	private PrintStream outstream=System.err;
	private static boolean verbose=false;
}
