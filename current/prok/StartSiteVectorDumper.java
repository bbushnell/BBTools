package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
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
 * Generates start-site training vectors for a dedicated StartNet.
 * For each NCBI CDS, outputs one-hot bases around the correct start (label=1)
 * and around each alternative start candidate from breakOrf (label=0).
 * Pure sequence context — no k-mer scores or heuristic features.
 * @author Brian Bushnell, Chloe
 */
public class StartSiteVectorDumper extends ProkObject {

	public static void main(String[] args){
		Timer t=new Timer();
		StartSiteVectorDumper x=new StartSiteVectorDumper(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public StartSiteVectorDumper(String[] args){
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
			}else if(a.equalsIgnoreCase("upstream") || a.equalsIgnoreCase("left")){
				upstream=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("downstream") || a.equalsIgnoreCase("right")){
				downstream=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("genomedir") || a.equalsIgnoreCase("dir")){
				assert(b!=null);
				for(String dir : b.split(",")){addGenomeDir(dir);}
			}else if(a.equalsIgnoreCase("scalars") || a.equalsIgnoreCase("addscalars")){
				includeScalars=parse.Parse.parseBoolean(b);
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
		assert(!fnaList.isEmpty()) : "At least 1 fasta file is required (in= or genomedir=).";
		assert(fnaList.size()==gffList.size()) : "Mismatched fna/gff counts.";
		assert(outFile!=null) : "Output file required (out=).";

		windowTotal=upstream+downstream;
		numScalarDims=(includeScalars ? 3 : 0);
		numInputDims=windowTotal*4+numScalarDims;
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

		ByteStreamWriter bsw=new ByteStreamWriter(outFile, true, false, false);
		bsw.start();

		ByteBuilder header=new ByteBuilder();
		header.append("#dims\t").append(numInputDims).append('\t').append(1).nl();
		bsw.print(header);

		long totalSites=0, totalPositive=0, totalNegative=0, totalGenes=0;

		for(int fnum=0; fnum<fnaList.size(); fnum++){
			String fna=fnaList.get(fnum);
			String gff=gffList.get(fnum);

			ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FASTA, -1);
			float gc=calcGC(reads);
			HashMap<String, HashSet<String>> ncbiStopToStarts=loadNcbiStartsByStop(gff);

			GeneCaller caller=new GeneCaller(minLen, 80, 110, -9999f, -9999f, -9999f, -9999f, -9999f, pgm);

			for(Read r : reads){
				ArrayList<Orf> orfs=caller.getAllScoredOrfs(r, pgm);
				byte[] bases=r.bases;
				String seqid=r.id;
				{
					int idx=seqid.indexOf(' ');
					if(idx>=0){seqid=seqid.substring(0, idx);}
				}

				HashMap<String, ArrayList<Orf>> byStop=new HashMap<String, ArrayList<Orf>>();
				for(Orf orf : orfs){
					if(orf.type!=CDS){continue;}
					String stopKey=seqid+":"+(orf.stop+1)+":"+orf.strand;
					ArrayList<Orf> group=byStop.get(stopKey);
					if(group==null){
						group=new ArrayList<Orf>();
						byStop.put(stopKey, group);
					}
					group.add(orf);
				}

				for(String stopKey : byStop.keySet()){
					HashSet<String> ncbiStarts=ncbiStopToStarts.get(stopKey);
					if(ncbiStarts==null){continue;}

					ArrayList<Orf> group=byStop.get(stopKey);
					totalGenes++;
					for(Orf orf : group){
						String startKey=String.valueOf(orf.start+1);
						boolean match=ncbiStarts.contains(startKey);

						ByteBuilder bb=new ByteBuilder();
						appendStartVector(bb, orf, bases, gc, match);
						bsw.print(bb);

						totalSites++;
						if(match){totalPositive++;}else{totalNegative++;}
					}
				}
			}
			if((fnum+1)%100==0){
				outstream.println("Processed "+(fnum+1)+"/"+fnaList.size()+" files, "+totalSites+" sites ("+totalPositive+" pos, "+totalNegative+" neg)");
			}
		}

		bsw.poisonAndWait();
		t.stop();
		outstream.println("Total genes: "+totalGenes);
		outstream.println("Total start sites: "+totalSites+" ("+totalPositive+" positive, "+totalNegative+" negative)");
		outstream.println("Positive rate: "+String.format("%.4f", totalPositive/(double)Tools.max(1, totalSites)));
		outstream.println(Tools.timeReadsBasesProcessed(t, fnaList.size(), totalSites, 8));
	}

	/** Load NCBI CDS annotations grouped by stop position.
	 * Returns map: "seqid:stop1based:strand" → set of "start1based" values. */
	private HashMap<String, HashSet<String>> loadNcbiStartsByStop(String gffFile){
		FileFormat ff=FileFormat.testInput(gffFile, FileFormat.GFF, null, true, true);
		ArrayList<GffLine>[] allLines=GffLine.loadGffFileByType(ff, "CDS", true);
		ArrayList<GffLine> cdsList=allLines[0];
		HashMap<String, HashSet<String>> map=new HashMap<String, HashSet<String>>();
		for(GffLine gl : cdsList){
			String seqid=gl.seqid;
			int idx=seqid.indexOf(' ');
			if(idx>=0){seqid=seqid.substring(0, idx);}
			String stopKey=seqid+":"+gl.stop+":"+gl.strand;
			HashSet<String> starts=map.get(stopKey);
			if(starts==null){
				starts=new HashSet<String>();
				map.put(stopKey, starts);
			}
			starts.add(String.valueOf(gl.start));
		}
		return map;
	}

	private void appendStartVector(ByteBuilder bb, Orf orf, byte[] bases, float gc, boolean match){
		if(orf.strand==Shared.PLUS){
			for(int offset=-upstream; offset<downstream; offset++){
				int pos=orf.start+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				appendOneHotBase(bb, x, offset>-upstream);
			}
		}else{
			for(int offset=upstream-1; offset>=-downstream; offset--){
				int pos=orf.stop+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				if(x>=0){x=3-x;}
				appendOneHotBase(bb, x, offset<upstream-1);
			}
		}
		if(includeScalars){
			bb.tab().append(orf.startScore, 6);
			bb.tab().append(gc, 4);
			bb.tab().append((float)(0.1*Math.log(Tools.max(1, orf.length()))), 6);
		}
		bb.tab().append(match ? 1 : 0);
		bb.nl();
	}

	private void appendOneHotBase(ByteBuilder bb, int x, boolean tab){
		if(tab){bb.tab();}
		bb.append(x==0 ? 1 : 0);
		bb.tab().append(x==1 ? 1 : 0);
		bb.tab().append(x==2 ? 1 : 0);
		bb.tab().append(x==3 ? 1 : 0);
	}

	private static float calcGC(ArrayList<Read> reads){
		long gc=0, at=0;
		for(Read r : reads){
			for(byte b : r.bases){
				int x=AminoAcid.baseToNumber[b];
				if(x==1 || x==2){gc++;}
				else if(x==0 || x==3){at++;}
			}
		}
		return (float)(gc/(double)Tools.max(1, gc+at));
	}

	/*--------------------------------------------------------------*/

	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> gffList=new ArrayList<String>();
	private String pgmFile=null;
	private String outFile=null;

	private int minLen=80;
	private int upstream=21;
	private int downstream=9;
	private boolean includeScalars=true;

	private final int windowTotal;
	private final int numScalarDims;
	private final int numInputDims;

	private PrintStream outstream=System.err;
}
