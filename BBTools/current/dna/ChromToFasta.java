package dna;

import java.io.File;
import java.util.ArrayList;

import fileIO.TextStreamWriter;
import shared.PreParser;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Jul 26, 2012
 *
 */
public class ChromToFasta {
	
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		if(args[0].contains("=") && (args[0].startsWith("build") || args[0].startsWith("genome"))){
			int build=Integer.parseInt(args[0].split("=")[1]);
//			Data.setGenome(build);
//			String s="", comma="";
//			for(int i=1; i<Data.numChroms; i++){
//				s=s+comma+Data.chromFname(i, build);
//				comma=",";
//			}
//			args[0]=s;
			args[0]=Data.chromFname(1, build);
			args[0]=args[0].substring(0, args[0].lastIndexOf('/'));
		}
		
		String[] chromfiles=args[0].split(",");
		
		if(chromfiles.length==1){
			File f=new File(chromfiles[0]);
			if(f.isDirectory()){
				ArrayList<String> list=new ArrayList<String>(4);
				for(File f2 : f.listFiles()){
					if(!f2.isDirectory() && f2.isFile()){
						String s=f2.getAbsolutePath();
						if(s.endsWith(".chrom") || s.endsWith(".chromC") || s.contains(".chrom.") || s.contains(".chromC.")){
							list.add(s);
						}
					}
				}
				chromfiles=list.toArray(new String[list.size()]);
			}
		}
		
		String outfile=args[1];
		int blocklen=Integer.parseInt(args[2]);
		int trigger=(args.length>3 ? Integer.parseInt(args[3]) : 0);
		
		TextStreamWriter tsw=new TextStreamWriter(outfile, true, false, false);
		tsw.start();
		
		if(trigger<=0){ //Write normally
			for(int i=0; i<chromfiles.length; i++){
				ChromosomeArray cha=ChromosomeArray.read(chromfiles[i]);
				writeChrom(cha, tsw, blocklen);
			}
		}else{ //Break into contigs
			int contig=1;
			for(int i=0; i<chromfiles.length; i++){
				ChromosomeArray cha=ChromosomeArray.read(chromfiles[i]);
				contig=writeContigs(cha, contig, trigger, blocklen, tsw);
			}
		}
		
		tsw.poison();
		
		try {tsw.join();}
		catch (InterruptedException e) {e.printStackTrace();}
		
		t.stop();
		System.err.println("Time:\t"+t);
	}
	
	public static int writeContigs(ChromosomeArray cha, int contig, int trigger, int fastaBlocklen, TextStreamWriter tsw){
		
		StringBuilder sb=new StringBuilder(4000);
		
		int ns=0;
		
		for(int aloc=cha.minIndex; aloc<=cha.maxIndex; aloc++){
			byte b=cha.get(aloc);
			if(b=='N'){
				ns++;
				if(sb.length()>0){
					sb.append('N');
					if(ns==trigger){
						sb.setLength(sb.length()-ns);
						tsw.print(">"+contig+"\n");
						contig++;
						writeContig(sb, tsw, fastaBlocklen);
						sb.setLength(0);
					}
				}
			}else{
				sb.append((char)b);
				ns=0;
			}
		}
		
		
		if(sb.length()>0){
			sb.setLength(sb.length()-ns);
			tsw.print(">"+contig+"\n");
			contig++;
			writeContig(sb, tsw, fastaBlocklen);
			sb.setLength(0);
		}
		
		return contig;
	}
	
	public static void writeContig(StringBuilder sb, TextStreamWriter tsw, int blocklen){
		for(int i=0; i<sb.length(); i+=blocklen){
			int max=Tools.min(i+blocklen, sb.length());
			tsw.println(sb.substring(i, max));
		}
	}
	
	public static void writeChrom(ChromosomeArray cha, String fname, int blocklen){
		TextStreamWriter tsw=new TextStreamWriter(fname, true, false, false);
		tsw.start();
		tsw.print(">"+cha.chromosome+"\n");
		writeChrom(cha, tsw, blocklen);
		tsw.poison();
	}
	
	public static void writeChrom(ChromosomeArray cha, TextStreamWriter tsw, int blocklen){
		tsw.println(">"+cha.chromosome);
		for(int i=0; i<=cha.maxIndex; i+=blocklen){
			int max=Tools.min(i+blocklen-1, cha.maxIndex);
			tsw.println(cha.getString(i, max));
		}
	}
	
}
