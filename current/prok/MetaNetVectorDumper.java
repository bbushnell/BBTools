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
import ml.CellNet;
import ml.CellNetParser;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/**
 * Generates training vectors for the MetaNet — the hierarchical network
 * that combines StartNet output, StopNet output, heuristic sub-scores,
 * codon bases, length, GC, and leader distance into a final ORF score adjustment.
 *
 * Input dims (34): StartNet(1), StopNet(1), heuristic startScore/stopScore/kmerScore/avgKmerScore(4),
 * log10(length)(1), GC(1), start codon one-hot bases(12), stop codon one-hot bases(12),
 * log10(leader distance)(1), orfLength/maxLength(1).
 *
 * @author Brian Bushnell, Chloe
 */
public class MetaNetVectorDumper extends ProkObject {

	public static void main(String[] args){
		Timer t=new Timer();
		MetaNetVectorDumper x=new MetaNetVectorDumper(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public MetaNetVectorDumper(String[] args){
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
			}else if(a.equalsIgnoreCase("startnets")){
				startNetPaths=b.split(",");
			}else if(a.equalsIgnoreCase("stopnets")){
				stopNetPaths=b.split(",");
			}else if(a.equalsIgnoreCase("gcmeans")){
				String[] parts=b.split(",");
				gcMeans=new float[parts.length];
				for(int j=0; j<parts.length; j++){gcMeans[j]=Float.parseFloat(parts[j]);}
			}else if(a.equalsIgnoreCase("genomedir") || a.equalsIgnoreCase("dir")){
				assert(b!=null);
				for(String dir : b.split(",")){addGenomeDir(dir);}
			}else if(a.equalsIgnoreCase("minLen") || a.equals("minlength")){
				minLen=Integer.parseInt(b);
			}else if(ProkObject.parse(arg, a, b)){
			}else if(parser.parse(arg, a, b)){
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if(pgmFile==null){pgmFile=Data.findPath("?model.pgm");}
		assert(!fnaList.isEmpty()) : "At least 1 fasta file required.";
		assert(fnaList.size()==gffList.size()) : "Mismatched fna/gff counts.";
		assert(outFile!=null) : "Output file required (out=).";
		assert(startNetPaths!=null) : "StartNet paths required (startnets=).";
		assert(stopNetPaths!=null) : "StopNet paths required (stopnets=).";
		assert(gcMeans!=null) : "GC means required (gcmeans=).";
		assert(startNetPaths.length==gcMeans.length) : "StartNet count must match gcmeans count.";
		assert(stopNetPaths.length==gcMeans.length) : "StopNet count must match gcmeans count.";

		startNets=new CellNet[startNetPaths.length];
		stopNets=new CellNet[stopNetPaths.length];
		for(int j=0; j<startNetPaths.length; j++){
			startNets[j]=CellNetParser.load(startNetPaths[j]);
			stopNets[j]=CellNetParser.load(stopNetPaths[j]);
		}
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

	private int selectNet(float gc){
		int best=0;
		float bestDist=Math.abs(gc-gcMeans[0]);
		for(int i=1; i<gcMeans.length; i++){
			float dist=Math.abs(gc-gcMeans[i]);
			if(dist<bestDist){bestDist=dist; best=i;}
		}
		return best;
	}

	void process(Timer t){
		ArrayList<String> pgmList=new ArrayList<String>();
		pgmList.add(pgmFile);
		final GeneModel pgm=PGMTools.loadAndMerge(pgmList);

		ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		ByteBuilder header=new ByteBuilder();
		header.append("#dims\t").append(NUM_INPUTS).append('\t').append(1).nl();
		bsw.print(header);

		long totalOrfs=0, totalPositive=0, totalNegative=0;

		for(int fnum=0; fnum<fnaList.size(); fnum++){
			String fna=fnaList.get(fnum);
			String gff=gffList.get(fnum);

			ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FASTA, -1);
			HashSet<String> ncbiCds=loadNcbiCds(gff);
			GeneCaller caller=new GeneCaller(minLen, 80, 110, -9999f, -9999f, -9999f, -9999f, -9999f, pgm);

			for(Read r : reads){
				byte[] bases=r.bases;
				float scaffoldGC=calcScaffoldGC(bases);
				int netIdx=selectNet(scaffoldGC);
				CellNet startNet=startNets[netIdx];
				CellNet stopNet=stopNets[netIdx];

				ArrayList<Orf> orfs=caller.getAllScoredOrfs(r, pgm);
				int[] maxLenByStop=computeMaxLenByStop(orfs);

				for(Orf orf : orfs){
					if(orf.type!=CDS){continue;}
					String key=makeKey(r.id, orf);
					boolean match=ncbiCds.contains(key);

					int maxLen=maxLenByStop[orf.stop];
					float startNetOut=runStartNet(startNet, orf, bases, scaffoldGC);
					float stopNetOut=runStopNet(stopNet, orf, bases, scaffoldGC);

					ByteBuilder bb=new ByteBuilder();
					appendMetaVector(bb, orf, startNetOut, stopNetOut, scaffoldGC, maxLen, match);
					bsw.print(bb);

					totalOrfs++;
					if(match){totalPositive++;}else{totalNegative++;}
				}
				scaffoldsProcessed++;
			}
			if((fnum+1)%100==0){
				outstream.println("Processed "+(fnum+1)+"/"+fnaList.size()+" files, "+scaffoldsProcessed+" scaffolds, "+totalOrfs+" orfs ("+totalPositive+" pos, "+totalNegative+" neg)");
			}
		}

		bsw.poisonAndWait();
		t.stop();
		outstream.println("Total ORFs: "+totalOrfs+" ("+totalPositive+" positive, "+totalNegative+" negative)");
		outstream.println("Positive rate: "+String.format("%.4f", totalPositive/(double)Tools.max(1, totalOrfs)));
	}

	private float runStartNet(CellNet net, Orf orf, byte[] bases, float gc){
		float[] in=new float[net.numInputs()];
		int idx=0;
		if(orf.strand==Shared.PLUS){
			for(int offset=-21; offset<9; offset++){
				int pos=orf.start+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
			}
		}else{
			for(int offset=20; offset>=-9; offset--){
				int pos=orf.stop+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				if(x>=0){x=3-x;}
				in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
			}
		}
		if(idx<in.length){
			in[idx++]=orf.startScore;
			in[idx++]=gc;
			in[idx++]=(float)(0.1*Math.log(Tools.max(1, orf.length())));
		}
		net.applyInput(in);
		return net.feedForward();
	}

	private float runStopNet(CellNet net, Orf orf, byte[] bases, float gc){
		float[] in=new float[net.numInputs()];
		int idx=0;
		int stopStart=orf.stop-2;
		if(orf.strand==Shared.PLUS){
			for(int offset=-9; offset<21; offset++){
				int pos=stopStart+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
			}
		}else{
			int stopAnchor=orf.start+2;
			for(int offset=8; offset>=-21; offset--){
				int pos=stopAnchor+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				if(x>=0){x=3-x;}
				in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
			}
		}
		if(idx<in.length){
			in[idx++]=orf.stopScore;
			in[idx++]=gc;
			in[idx++]=(float)(0.1*Math.log(Tools.max(1, orf.length())));
		}
		net.applyInput(in);
		return net.feedForward();
	}

	private void appendMetaVector(ByteBuilder bb, Orf orf, float startNetOut, float stopNetOut,
			float gc, int maxLen, boolean match){
		bb.append(startNetOut, 6);
		bb.tab().append(stopNetOut, 6);
		bb.tab().append(orf.startScore, 6);
		bb.tab().append(orf.stopScore, 6);
		bb.tab().append(orf.kmerScore, 6);
		bb.tab().append(orf.averageKmerScore(), 6);
		bb.tab().append((float)Math.log10(Tools.max(1, orf.length())), 6);
		bb.tab().append(gc, 4);

		appendCodonOneHot(bb, orf.startCodon);
		appendCodonOneHot(bb, orf.stopCodon);

		float leaderDist=(maxLen>0 ? maxLen : orf.length());
		bb.tab().append((float)Math.log10(Tools.max(1, leaderDist)), 6);
		bb.tab().append(orf.length()/(float)Tools.max(1, leaderDist), 6);

		bb.tab().append(match ? 1 : 0);
		bb.nl();
	}

	private void appendCodonOneHot(ByteBuilder bb, int packedCodon){
		for(int pos=2; pos>=0; pos--){
			int base=(packedCodon>>(pos*2))&3;
			bb.tab().append(base==0?1:0);
			bb.tab().append(base==1?1:0);
			bb.tab().append(base==2?1:0);
			bb.tab().append(base==3?1:0);
		}
	}

	private int[] computeMaxLenByStop(ArrayList<Orf> orfs){
		int maxStop=0;
		for(Orf orf : orfs){if(orf.stop>maxStop){maxStop=orf.stop;}}
		int[] maxLen=new int[maxStop+1];
		for(Orf orf : orfs){
			if(orf.type==CDS && orf.length()>maxLen[orf.stop]){
				maxLen[orf.stop]=orf.length();
			}
		}
		return maxLen;
	}

	private HashSet<String> loadNcbiCds(String gffFile){
		FileFormat ff=FileFormat.testInput(gffFile, FileFormat.GFF, null, true, true);
		ArrayList<GffLine>[] allLines=GffLine.loadGffFileByType(ff, "CDS", true);
		ArrayList<GffLine> cdsList=allLines[0];
		HashSet<String> set=new HashSet<String>(cdsList.size()*2);
		for(GffLine gl : cdsList){
			String seqid=gl.seqid;
			int idx=seqid.indexOf(' ');
			if(idx>=0){seqid=seqid.substring(0, idx);}
			set.add(seqid+":"+gl.start+":"+gl.stop+":"+gl.strand);
		}
		return set;
	}

	private String makeKey(String readId, Orf orf){
		String seqid=readId;
		int idx=seqid.indexOf(' ');
		if(idx>=0){seqid=seqid.substring(0, idx);}
		return seqid+":"+(orf.start+1)+":"+(orf.stop+1)+":"+orf.strand;
	}

	private static float calcScaffoldGC(byte[] bases){
		long gc=0, at=0;
		for(byte b : bases){
			int x=AminoAcid.baseToNumber[b];
			if(x==1 || x==2){gc++;}
			else if(x==0 || x==3){at++;}
		}
		return (float)(gc/(double)Tools.max(1, gc+at));
	}

	/*--------------------------------------------------------------*/

	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> gffList=new ArrayList<String>();
	private String pgmFile=null;
	private String outFile=null;
	private String[] startNetPaths=null;
	private String[] stopNetPaths=null;
	private float[] gcMeans=null;
	private CellNet[] startNets;
	private CellNet[] stopNets;
	private int minLen=80;
	private long scaffoldsProcessed=0;

	static final int NUM_INPUTS=34;

	private PrintStream outstream=System.err;
}
