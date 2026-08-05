package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import dna.AminoAcid;
import dna.Data;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import gff.GffLine;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/**
 * Generates training vectors for the CallGenes neural network.
 * For each ORF candidate on each scaffold, writes a feature vector with
 * an NCBI-match label. Output is tab-delimited, compatible with ml.Trainer.
 * @author Brian Bushnell, Chloe
 */
public class GeneVectorDumper extends ProkObject {

	public static void main(String[] args){
		Timer t=new Timer();
		GeneVectorDumper x=new GeneVectorDumper(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public GeneVectorDumper(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(PGMTools.parseStatic(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("infna") || a.equals("fna") || a.equals("ref")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(a.equalsIgnoreCase("ingff") || a.equalsIgnoreCase("gff") || a.equalsIgnoreCase("gffin")){
				assert(b!=null);
				Tools.addFiles(b, gffList);
			}else if(a.equals("pgm") || a.equals("gm") || a.equals("model")){
				if(b!=null){
					if(b.equalsIgnoreCase("auto") || b.equalsIgnoreCase("default")){
						b=Data.findPath("?model.pgm");
					}
					pgmFile=b;
				}
			}else if(a.equals("out")){
				outFile=b;
			}else if(a.equalsIgnoreCase("startWindow") || a.equalsIgnoreCase("sw")){
				startWindowLeft=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("startWindowRight") || a.equalsIgnoreCase("swr")){
				startWindowRight=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("stopWindow") || a.equalsIgnoreCase("stw")){
				stopWindowLeft=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("stopWindowRight") || a.equalsIgnoreCase("stwr")){
				stopWindowRight=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("onehot") || a.equalsIgnoreCase("bases")){
				includeOneHot=parse.Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("genomedir") || a.equalsIgnoreCase("dir")){
				assert(b!=null);
				for(String dir : b.split(",")){
					addGenomeDir(dir);
				}
			}else if(a.equalsIgnoreCase("minLen") || a.equals("minlength")){
				minLen=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("startSlop")){
				startSlop=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("stopSlop")){
				stopSlop=Integer.parseInt(b);
			}else if(ProkObject.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if(pgmFile==null){pgmFile=Data.findPath("?model.pgm");}
		assert(!fnaList.isEmpty()) : "At least 1 fasta file is required (in= or genomedir=).";
		assert(fnaList.size()==gffList.size()) : "Must have equal numbers of fna and gff files; got "+fnaList.size()+" fna, "+gffList.size()+" gff.";
		assert(outFile!=null) : "Output file required (out=).";

		startWindowTotal=startWindowLeft+startWindowRight;
		stopWindowTotal=stopWindowLeft+stopWindowRight;
		numOneHotDims=(startWindowTotal+stopWindowTotal)*4;
		numScalarDims=15;
		numInputDims=numScalarDims+(includeOneHot ? numOneHotDims : 0);
	}

	private void addGenomeDir(String dirPath){
		File dir=new File(dirPath);
		assert(dir.isDirectory()) : "Not a directory: "+dirPath;
		File[] files=dir.listFiles();
		Arrays.sort(files);
		int added=0;
		for(File f : files){
			String name=f.getName();
			if(name.endsWith(".fna.gz") || name.endsWith(".fna")){
				String gffName=name.replace(".fna.gz", ".gff.gz").replace(".fna", ".gff");
				File gff=new File(dir, gffName);
				if(gff.exists()){
					fnaList.add(f.getAbsolutePath());
					gffList.add(gff.getAbsolutePath());
					added++;
				}
			}
		}
		outstream.println("Added "+added+" genome/gff pairs from "+dirPath);
	}

	void process(Timer t){
		ArrayList<String> pgmList=new ArrayList<String>();
		pgmList.add(pgmFile);
		final GeneModel pgm=PGMTools.loadAndMerge(pgmList);
		final float gc=pgm.gc();

		ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		ByteBuilder header=new ByteBuilder();
		header.append("#dims\t").append(numInputDims).append('\t').append(1).append('\t').append(1).nl();
		header.append("#features\tstartScore\tstopScore\tkmerScore\tavgKmerScore\tlog10len");
		header.append("\tgc\tframe0\tframe1\tframe2\tstartATG\tstartGTG\tstartTTG");
		header.append("\tstopTAG\tstopTAA\tstopTGA");
		if(includeOneHot){
			for(int i=0; i<startWindowTotal; i++){
				header.append("\tSb").append(i).append("A\tSb").append(i).append("C\tSb").append(i).append("G\tSb").append(i).append("T");
			}
			for(int i=0; i<stopWindowTotal; i++){
				header.append("\tEb").append(i).append("A\tEb").append(i).append("C\tEb").append(i).append("G\tEb").append(i).append("T");
			}
		}
		header.append("\tweight\ttarget").nl();
		bsw.print(header);

		long totalOrfs=0, totalPositive=0, totalNegative=0;

		for(int fnum=0; fnum<fnaList.size(); fnum++){
			String fna=fnaList.get(fnum);
			String gff=gffList.get(fnum);

			ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FASTA, -1);
			HashSet<String> ncbiCds=loadNcbiCds(gff);

			GeneCaller caller=new GeneCaller(minLen, 80, 110, -0.10f, -0.5f, 0.02f, -9999f, -9999f, pgm);

			for(Read r : reads){
				ArrayList<Orf> orfs=caller.getAllScoredOrfs(r, pgm);
				byte[] bases=r.bases;

				for(Orf orf : orfs){
					String key=makeKey(r.id, orf);
					boolean match=ncbiCds.contains(key);

					ByteBuilder bb=new ByteBuilder();
					appendVector(bb, orf, bases, gc, match);
					bsw.print(bb);

					totalOrfs++;
					if(match){totalPositive++;}else{totalNegative++;}
				}
			}
			if((fnum+1)%100==0){
				outstream.println("Processed "+(fnum+1)+"/"+fnaList.size()+" files, "+totalOrfs+" orfs ("+totalPositive+" pos, "+totalNegative+" neg)");
			}
		}

		bsw.poisonAndWait();
		t.stop();
		outstream.println("Total ORFs: "+totalOrfs+" ("+totalPositive+" positive, "+totalNegative+" negative)");
		outstream.println("Positive rate: "+String.format("%.4f", totalPositive/(double)Tools.max(1, totalOrfs)));
		outstream.println(Tools.timeReadsBasesProcessed(t, fnaList.size(), totalOrfs, 8));
	}

	/** Load NCBI CDS annotations into a set of "seqid:start:stop:strand" keys. */
	private HashSet<String> loadNcbiCds(String gffFile){
		FileFormat ff=FileFormat.testInput(gffFile, FileFormat.GFF, null, true, true);
		ArrayList<GffLine>[] allLines=GffLine.loadGffFileByType(ff, "CDS", true);
		ArrayList<GffLine> cdsList=allLines[0];
		HashSet<String> set=new HashSet<String>(cdsList.size()*2);
		for(GffLine gl : cdsList){
			String seqid=gl.seqid;
			int idx=seqid.indexOf(' ');
			if(idx>=0){seqid=seqid.substring(0, idx);}
			for(int sOff=-startSlop; sOff<=startSlop; sOff++){
				for(int eOff=-stopSlop; eOff<=stopSlop; eOff++){
					set.add(seqid+":"+(gl.start+sOff)+":"+(gl.stop+eOff)+":"+gl.strand);
				}
			}
		}
		return set;
	}

	/** Create a matching key for an Orf (converting 0-based to 1-based). */
	private String makeKey(String readId, Orf orf){
		String seqid=readId;
		int idx=seqid.indexOf(' ');
		if(idx>=0){seqid=seqid.substring(0, idx);}
		return seqid+":"+(orf.start+1)+":"+(orf.stop+1)+":"+orf.strand;
	}

	private void appendVector(ByteBuilder bb, Orf orf, byte[] bases, float gc, boolean match){
		float startScore=orf.startScore;
		float stopScore=orf.stopScore;
		float kmerScore=orf.kmerScore;
		float avgKmerScore=orf.averageKmerScore();
		float logLen=(float)Math.log10(Tools.max(1, orf.length()));
		int frame=orf.frame;

		int startCodonType=classifyStartCodon(orf.startCodon);
		int stopCodonType=classifyStopCodon(orf.stopCodon);

		bb.append(startScore, 6);
		bb.tab().append(stopScore, 6);
		bb.tab().append(kmerScore, 6);
		bb.tab().append(avgKmerScore, 6);
		bb.tab().append(logLen, 6);
		bb.tab().append(gc, 4);

		for(int f=0; f<3; f++){bb.tab().append(f==frame ? 1 : 0);}
		for(int c=0; c<3; c++){bb.tab().append(c==startCodonType ? 1 : 0);}
		for(int c=0; c<3; c++){bb.tab().append(c==stopCodonType ? 1 : 0);}

		if(includeOneHot){
			if(orf.strand==Shared.PLUS){
				appendOneHotBases(bb, bases, orf.start, -startWindowLeft, startWindowRight, false);
				appendOneHotBases(bb, bases, orf.stop-2, -stopWindowLeft, stopWindowRight, false);
			}else{
				appendOneHotBases(bb, bases, orf.stop, -startWindowRight+1, startWindowLeft+1, true);
				appendOneHotBases(bb, bases, orf.start+2, -stopWindowRight+1, stopWindowLeft+1, true);
			}
		}

		bb.tab().append(1.0f, 1);
		bb.tab().append(match ? 1 : 0);
		bb.nl();
	}

	private void appendOneHotBases(ByteBuilder bb, byte[] bases, int anchor, int relLeft, int relRight, boolean reverseComplement){
		if(!reverseComplement){
			for(int offset=relLeft; offset<relRight; offset++){
				int pos=anchor+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				appendOneHotBase(bb, x);
			}
		}else{
			for(int offset=relRight-1; offset>=relLeft; offset--){
				int pos=anchor+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				if(x>=0){x=3-x;}
				appendOneHotBase(bb, x);
			}
		}
	}

	private void appendOneHotBase(ByteBuilder bb, int x){
		bb.tab().append(x==0 ? 1 : 0);
		bb.tab().append(x==1 ? 1 : 0);
		bb.tab().append(x==2 ? 1 : 0);
		bb.tab().append(x==3 ? 1 : 0);
	}

	/** ATG=0, GTG=1, TTG=2 */
	private static int classifyStartCodon(int packed){
		String s=AminoAcid.codonToString(packed);
		if("ATG".equals(s)){return 0;}
		if("GTG".equals(s)){return 1;}
		if("TTG".equals(s)){return 2;}
		return 0;
	}

	/** TAG=0, TAA=1, TGA=2 */
	private static int classifyStopCodon(int packed){
		String s=AminoAcid.codonToString(packed);
		if("TAG".equals(s)){return 0;}
		if("TAA".equals(s)){return 1;}
		if("TGA".equals(s)){return 2;}
		return 0;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> gffList=new ArrayList<String>();
	private String pgmFile=null;
	private String outFile=null;

	private int minLen=80;
	private int startSlop=0;
	private int stopSlop=0;

	private int startWindowLeft=20;
	private int startWindowRight=10;
	private int stopWindowLeft=10;
	private int stopWindowRight=5;
	private boolean includeOneHot=false;

	private final int startWindowTotal;
	private final int stopWindowTotal;
	private final int numOneHotDims;
	private final int numScalarDims;
	private final int numInputDims;

	private PrintStream outstream=System.err;
}
