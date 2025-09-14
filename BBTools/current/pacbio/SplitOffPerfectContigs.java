package pacbio;

import java.io.File;
import java.util.ArrayList;

import dna.ChromosomeArray;
import dna.Data;
import dna.FastaToChromArrays2;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Timer;
import shared.Tools;
import stream.SiteScore;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.Range;

/**
 * @author Brian Bushnell
 * @date Jul 26, 2012
 *
 */
public class SplitOffPerfectContigs {
	
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
//		ChromosomeArray c=new ChromosomeArray(1, (byte)1, "ANNNAAAANAAANNA");
//		System.out.println(c.toContigRanges(3));
//		System.out.println(c.toContigRanges(2));
//		System.out.println(c.toContigRanges(1));
//		assert(false);
		
		
		Data.GENOME_BUILD=-1;
		String dest=null;
		String covfile=null;
		String sitesfile=null;
		String contigfile=null;
		int trigger=50;
		int blocklen=100;
		int mincoverage=2;
		int padding=4;
		int buildout=-1;
		String name=null;
		String source=null;
		
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("genome") || a.equals("build")){
				Data.setGenome(Integer.parseInt(b));
				name=Data.name;
				source=Data.genomeSource;
				System.out.println("Set Data.GENOME_BUILD to "+Data.GENOME_BUILD);
			}else if(a.equals("outgenome") || a.equals("outbuild") || a.equals("genomeout") || a.equals("buildout")){
				buildout=Integer.parseInt(b);
			}else if(a.equals("out") || a.equals("outfile")){
				dest=b;
			}else if(a.startsWith("cov") || a.startsWith("pcov") || a.startsWith("perfectcov")){
				covfile=b;
			}else if(a.startsWith("sites") || a.startsWith("psites") || a.startsWith("perfectsites")){
				sitesfile=b;
			}else if(a.equals("padding")){
				padding=Integer.parseInt(b);
			}else if(a.equals("trigger")){
				trigger=Integer.parseInt(b);
			}else if(a.startsWith("mincov")){
				mincoverage=Integer.parseInt(b);
			}else if(a.equals("blocklen")){
				blocklen=Integer.parseInt(b);
			}else if(a.equals("contigfile")){
				contigfile=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.startsWith("breakbad") || a.startsWith("splitbad") || a.startsWith("splitchim")){
				BREAK_BAD_CONTIGS=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter: "+args[i]);
			}
		}
		
		assert(Data.GENOME_BUILD>-1);
		if(buildout<=0){buildout=Data.GENOME_BUILD;}
//		assert(buildout!=Data.GENOME_BUILD); //For testing
		
		TextStreamWriter tsw=new TextStreamWriter(dest, false, true, false);
		tsw.start();
		
		//Break into contigs
		long contig=1;
		
		if(contigfile!=null){
			if(new File(contigfile).exists()){
				TextFile tf=new TextFile(contigfile, false);
				String s=tf.nextLine();
				if(s!=null){contig=Long.parseLong(s);}
				tf.close();
			}
		}
		
		ArrayList<CoverageArray> calist=null;
		if(sitesfile!=null){
			calist=toCoverage(sitesfile, padding);
			System.out.println("Made coverage; list size is "+calist.size());
		}
		
		if(buildout==Data.GENOME_BUILD){
			String fname=Data.chromFname(1, buildout);
			fname=fname.replaceFirst("/genome/", "/index/");
			fname=fname.substring(0, fname.lastIndexOf('/'));
			File dir=new File(fname);
			if(dir.exists()){
				System.out.println("Deleting old index.");
				for(File f2 : dir.listFiles()){
					if(f2.isFile() && !f2.isDirectory() && f2.getName().contains(".int2d")){f2.delete();}
				}
			}
		}
		
		for(int chrom=1; chrom<=Data.numChroms; chrom++){
			ChromosomeArray cha=Data.getChromosome(chrom);
			Data.unload(chrom, true);
			CoverageArray ca=null;
			if(calist!=null){
				if(calist.size()>chrom){
					ca=calist.get(chrom);
					calist.set(chrom, null);
				}
			}else{
				assert(covfile!=null && covfile.contains("#"));
				ca=ReadWrite.read(CoverageArray.class, covfile.replaceFirst("#", ""+chrom), true);
				if(ca==null){System.out.println("Can't find coverage for chrom "+chrom+" in file "+covfile.replaceFirst("#", ""+chrom));}
			}
			if(ca!=null){
				contig=writeContigs(cha, ca, contig, trigger, mincoverage, blocklen, tsw, buildout, shared.Tools.max(1, padding));
			}else{
				System.out.println("Can't find coverage for chrom "+chrom);
			}
		}
		
		
		tsw.poison();
		
		if(contigfile!=null){
			ReadWrite.writeString(""+contig, contigfile, false);
		}
		
		FastaToChromArrays2.writeInfo(buildout, Data.numChroms, name, source, false, false);
		
		t.stop();
		
		System.out.println("          \tWrote   \tKept      \tDropped   \tSplit");
		System.out.println("Bases     \t"+basesWritten+"   \t"+basesKept+"   \t"+basesDropped+"       \t"+basesX);
		System.out.println("Contigs   \t"+contigsWritten+"     \t"+contigsKept+"     \t"+contigsDropped+"        \t"+contigsX);
		System.out.println("Avg Len   \t"+(basesWritten/Tools.max(contigsWritten,1))+"     \t"+(basesKept/Tools.max(contigsKept,1))
				+"     \t"+(basesDropped/Tools.max(contigsDropped, 1))+"        \t"+(basesX/Tools.max(contigsX, 1)));
		
		System.out.println("Time:\t"+t);
	}
	
	public static long writeContigs(ChromosomeArray cha, CoverageArray ca, long contig, int trigger, int minAcceptableCoverage, int fastaBlocklen,
			TextStreamWriter tsw, int buildout, int tipbuffer){
		
		ArrayList<Range> list=cha.toContigRanges(trigger);
		
		int minContig=MIN_CONTIG_TO_ADD;
		
		if(BREAK_BAD_CONTIGS){
			for(Range r : list){
				if(r.length()>=minContig){
//					int uncovered=0;
//					for(int i=r.a; i<=r.b; i++){
//						int cov=ca.get(i);
//						if(cov<minAcceptableCoverage){uncovered++;}
//					}
					
					
					//Forward pass
					int lastx=-1000;
					int contiglen=0;
					for(int i=r.a; i<=r.b; i++){
						int cov=ca.get(i);
						if(cov<minAcceptableCoverage){
							if(contiglen>=minContig){
								byte c=cha.get(i);
								if(c!='N' && c!='X'){basesX++;}
								if(i-lastx>10){
									contigsX++;
								}
								cha.set(i, 'X');
								lastx=i;
							}
							contiglen=0;
						}else{
							contiglen++;
						}
					}
					
					//Reverse pass
					lastx=Integer.MAX_VALUE;
					contiglen=0;
					for(int i=r.b; i>=r.a; i--){
						int cov=ca.get(i);
						if(cov<minAcceptableCoverage){
							if(contiglen>=minContig){
								byte c=cha.get(i);
								if(c!='N' && c!='X'){basesX++;}
								if(lastx-i>10){
									contigsX++;
								}
								cha.set(i, 'X');
								lastx=i;
							}
							contiglen=0;
						}else{
							contiglen++;
						}
					}
				}
			}
			list=cha.toContigRanges(trigger);
		}
		
		
		ArrayList<Range> good=new ArrayList<Range>();
		ArrayList<Range> bad=new ArrayList<Range>();
		int badlen=0;
		
		for(Range r : list){
			final int length=r.length();
			if(length>=minContig){
				int minCov=Integer.MAX_VALUE;
				for(int i=r.a+tipbuffer; i<=r.b-tipbuffer; i++){
					minCov=Tools.min(minCov, ca.get(i));
				}
				if(minCov>=minAcceptableCoverage){
					good.add(r);
					if(verbose){
						StringBuilder sb0=new StringBuilder(), sb1=new StringBuilder(), sb2=new StringBuilder();
						for(int i=r.a; i<=r.b; i++){
							int cov=ca.get(i);
							char b=(char) cha.get(i);
							sb0.append(b);
							sb1.append(b+"\t");
							sb2.append(cov+"\t");
						}
						System.out.println(sb0+"\n"+sb1+"\n"+sb2+"\n");
					}
				}else{
					bad.add(r);
					badlen+=length+N_PAD_LENGTH;
					if(verbose){
						StringBuilder sb0=new StringBuilder(), sb1=new StringBuilder(), sb2=new StringBuilder();
						for(int i=r.a; i<=r.b; i++){
							int cov=ca.get(i);
							char b=(char) cha.get(i);
							sb0.append(b);
							sb1.append(b+"\t");
							sb2.append(cov+"\t");
						}
						System.err.println(sb0+"\n"+sb1+"\n"+sb2+"\n");
					}
				}
			}else{
				contigsDropped++;
				basesDropped+=length;
			}
		}
		
		for(Range r : good){
			contigsWritten++;
			basesWritten+=r.length();
			String s=cha.getString(r.a, r.b);
			tsw.print(">"+contig+"\n");
			contig++;
			writeContig(s, tsw, fastaBlocklen);
//			for(int i=r.a; i<=r.b; i++){cha.set(i, 'N');} //Delete "good" contigs from reference.
		}
		
		badlen=badlen+2*N_PAD_LENGTH2-N_PAD_LENGTH+10;
		ChromosomeArray cha2=new ChromosomeArray(cha.chromosome, cha.strand, 0, badlen);
		cha2.maxIndex=-1;
		cha2.minIndex=0;
		for(int i=0; i<N_PAD_LENGTH2; i++){
			cha2.set(i, 'N');
		}
		for(Range r : bad){
			contigsKept++;
			basesKept+=r.length();
			
			String s=cha.getString(r.a, r.b);
			for(int i=0; i<s.length(); i++){
				cha2.set(cha2.maxIndex+1, s.charAt(i));
			}
			for(int i=0; i<N_PAD_LENGTH; i++){
				cha2.set(cha2.maxIndex+1, 'N');
			}
		}
		for(int i=N_PAD_LENGTH; i<N_PAD_LENGTH2; i++){
			cha2.set(cha2.maxIndex+1, 'N');
		}

//		ReadWrite.writeObjectInThread(cha2, Data.chromFname(cha2.chromosome, Data.GENOME_BUILD));
		String fname=Data.chromFname(cha2.chromosome, buildout);
		{
			File f=new File(fname.substring(0, fname.lastIndexOf('/')));
			if(!f.exists()){
				f.mkdirs();
			}
		}
		
		ReadWrite.write(cha2, fname, false);
		
		return contig;
	}
	
	public static void writeContig(CharSequence sb, TextStreamWriter tsw, int blocklen){
		for(int i=0; i<sb.length(); i+=blocklen){
			int max=Tools.min(i+blocklen, sb.length());
			tsw.println(sb.subSequence(i, max));
		}
	}
	
	
	public static ArrayList<CoverageArray> toCoverage(String sitesfile, int padding){
		ArrayList<CoverageArray> pcov=new ArrayList<CoverageArray>(8);
		pcov.add(new CoverageArray2(0,1000));
		
		long perfect=0;
		long semiperfect=0;
		long sites=0;
		
		String[] files=sitesfile.split(",");
		for(String f : files){
			TextFile tf=new TextFile(f, false);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				String[] split=line.split("\t");
				for(String s : split){
					SiteScore ss=SiteScore.fromText(s);
					while(pcov.size()<=ss.chrom){
						pcov.add(new CoverageArray2(pcov.size(), 500));
					}
					if(ss.perfect || ss.semiperfect){
						CoverageArray ca=pcov.get(ss.chrom);
						for(int i=ss.start+padding; i<=ss.stop-padding; i++){
							ca.increment(i);
						}
					}
					if(ss.perfect){perfect++;}
					if(ss.semiperfect){semiperfect++;}
					sites++;
					assert(!ss.perfect || ss.semiperfect) : ss.perfect+", "+ss.semiperfect+"\n"+ss.header()+"\n"+ss.toText()+"\n"+s+"\n";
				}
			}
			tf.close();
		}
		System.out.println("Read "+files.length+" sites file"+(files.length==1 ? "." : "s."));
		System.out.println("sites="+sites+"  \tsemiperfect="+semiperfect+"  \tperfect="+perfect);
		return pcov;
	}
	

	public static long basesWritten=0;
	public static long basesKept=0;
	public static long basesDropped=0;
	public static long basesX=0;
	public static long contigsWritten=0;
	public static long contigsKept=0;
	public static long contigsDropped=0;
	public static long contigsX=0;
	
	public static int N_PAD_LENGTH=MergeFastaContigs.N_PAD_LENGTH;
	public static int N_PAD_LENGTH2=MergeFastaContigs.N_PAD_LENGTH2; //for ends
	public static int MIN_CONTIG_TO_ADD=50;
	public static boolean BREAK_BAD_CONTIGS=false;
	
	public static boolean verbose=false;
	
}
