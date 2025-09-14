package jgi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Dec 13, 2012
 *
 */
public class CountGC {
	
	public static void main(String[] args){
		
		Timer t=new Timer();
		if(args.length==0){
			System.out.println("Usage: CountGC in=<infile> out=<outfile>");
			System.out.println("Alternately, 'out=stdout' will print to standard out.");
			System.out.println("Optional flag, format:");
			System.out.println("format=1\tid start stop A C G T N GC");
			System.out.println("format=2\tid gc");
			System.out.println("format=4\tid length gc");
			System.out.println("Output is always tab-delimited.  AGCT are fractions of defined bases; N is fraction of total bases.");
			System.exit(0);
		}
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean benchmark=false;
		ReadWrite.USE_UNPIGZ=true;
		
		String in=null, out=null;
		
		for(int i=0; i<args.length; i++){

			if(true){
				final String arg=args[i];
				final String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				
				if(Parser.parseCommonStatic(arg, a, b)){
					//do nothing
				}else if(Parser.parseZip(arg, a, b)){
					//do nothing
				}else if(Parser.parseQuality(arg, a, b)){
					//do nothing
				}else if(a.equals("in")){
					in=b;
				}else if(a.equals("out")){
					out=b;
					if(b==null || "summaryonly".equalsIgnoreCase(b) || "none".equalsIgnoreCase(b)){
						out=null;
						SUMMARY_ONLY=true;
					}else if("benchmark".equalsIgnoreCase(b)){
						benchmark=true;
						out=null;
						SUMMARY_ONLY=true;
					}
				}else if(a.equals("benchmark")){
					benchmark=Parse.parseBoolean(b);
					if(benchmark){
						out=null;
						SUMMARY_ONLY=true;
					}
				}else if(a.equals("format")){
					FORMAT=Integer.parseInt(b);
					if(FORMAT!=1 && FORMAT!=2 && FORMAT!=4){
						throw new RuntimeException("\nUnknown format: "+FORMAT+"; valid values are 1, 2, and 4.\n");
					}
				}else if(in==null && i==0 && !args[i].contains("=")){
					in=args[i];
				}else if(out==null && i==1 && !args[i].contains("=")){
					out=args[i];
				}
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		long[] counts=null;
		long sum=0;
		
		if(out==null || out.equalsIgnoreCase("stdout") || out.equalsIgnoreCase("standardout")){out=null;}
		
		InputStream is=null;
		{
			if(in==null){throw new RuntimeException("No input file.");}
			if(in.equalsIgnoreCase("stdin") || in.equalsIgnoreCase("standardin")){
				is=System.in;
			}else{
				File f=new File(in);
				if((!f.exists() || f.isDirectory()) && !in.toLowerCase().startsWith("stdin")){
					throw new RuntimeException("Input file does not appear to be valid: "+in);
				}
			}
		}
		
		if(is==null){is=ReadWrite.getInputStream(in, false, true);}
		try {
			if(benchmark){sum=bench2(is);}
			else{
				FileFormat ff=FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
				boolean fastq=ff.fastq();
				boolean fasta=!fastq; //Default.
				if(fastq){counts=countFastq(is, out);}
				else if(fasta){counts=countFasta(is, out);}
				else{throw new RuntimeException("Unknown or unsupported file format.");}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			if(is!=System.in){is.close();}
		} catch (IOException e) {
			e.printStackTrace();
		}
		

		t.stop();
		
		if(benchmark){
			outstream.println("Time: \t"+t);
			long bytes=new File(in).length();
			if(bytes<1){bytes=LIMSUM;}
			double mbps1=bytes*1000d/t.elapsed;
			double mbps2=sum*1000d/t.elapsed;
			outstream.println(Tools.format("Raw Speed:         \t%.2f MBytes/s",mbps1));
			outstream.println(Tools.format("Uncompressed Speed:\t%.2f MBytes/s",mbps2));
		}else{
			outstream.println(toString2(new StringBuilder("Overall"), counts));
			outstream.println("Time: \t"+t);
			long bytes=new File(in).length();
			if(bytes<1){bytes=LIMSUM;}
			double mbps=bytes*1000d/t.elapsed;
			double mbpps=shared.Vector.sum(counts)*1000d/t.elapsed;
			outstream.println(Tools.format("Speed:\t%.2f MBytes/s",mbps));
			outstream.println(Tools.format("      \t%.2f MBases/s",mbpps));
		}
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public static long bench2(InputStream is) throws IOException{
		final byte[] buf=new byte[32768];
		long sum=0;
		for(long len=is.read(buf); len>0; len=is.read(buf)){sum+=len;}
		return sum;
	}
	
	public static long[] countFasta(InputStream is, String out) throws IOException{
		
		long limsum=0;
		final byte[] buf=new byte[32768];
		final TextStreamWriter tsw=(out==null ? null : new TextStreamWriter(out, true, false, false));
		if(tsw!=null){tsw.start();}
		final int[] counts=KillSwitch.allocInt1D(6);
		final long[] overall=new long[6];
		final StringBuilder hdr=new StringBuilder();
		boolean hdmode=false;
		
		int i=0;
		int lim=is.read(buf);
		limsum+=lim;
		
		while(lim>0){
			if(hdmode){
				while(i<lim){
					byte c=buf[i];
					i++;
					if(c<=slashr){hdmode=false; break;}
					hdr.append((char)c);
				}
			}
			
			if(!hdmode){
				while(i<lim){
					byte c=buf[i];
					i++;
					
					if(c==carrot){
						hdmode=true;
						if(hdr.length()>0 || shared.Vector.sum(counts)>0){
							if(tsw!=null){tsw.print(toString2(hdr, counts));}else if(!SUMMARY_ONLY){System.out.print(toString2(hdr, counts));}
							hdr.setLength(0);
							for(int j=0; j<counts.length; j++){
								overall[j]+=counts[j];
								counts[j]=0;
							}
						}
						break;
					}
					counts[charToNum[c]]++;
					
				}
			}
			if(i>=lim){
				i=0;
				lim=is.read(buf);
				limsum+=lim;
			}
		}
		
		if(hdr.length()>0 || shared.Vector.sum(counts)>0){
			if(tsw!=null){tsw.print(toString2(hdr, counts));}else if(!SUMMARY_ONLY){System.out.println(toString2(hdr, counts));}
			hdr.setLength(0);
			for(int j=0; j<counts.length; j++){
				overall[j]+=counts[j];
				counts[j]=0;
			}
		}
		
		if(tsw!=null){
			tsw.poisonAndWait();
		}
		LIMSUM=limsum;
		return overall;
	}
	
	public static long[] countFastq(InputStream is, String out) throws IOException{
//		assert(false) : "Fastq mode - TODO"; //TODO
		long limsum=0;
		final byte[] buf=new byte[32768];
		final TextStreamWriter tsw=(out==null ? null : new TextStreamWriter(out, true, false, false));
		if(tsw!=null){tsw.start();}
		final int[] counts=KillSwitch.allocInt1D(6);
		final long[] overall=new long[6];
		final StringBuilder hdr=new StringBuilder();
		
		int mode=3;
		
		int i=0;
		int lim=is.read(buf);
		limsum+=lim;
		
		while(mode==3 && lim>0){
			while(i<lim && buf[i]!=at){i++;}
			if(i>=lim){
				lim=is.read(buf);
				limsum+=lim;
			}else{
				assert(buf[i]==at);
				mode=0;
			}
		}
		
		while(lim>0){
			if(mode==0){
				while(i<lim){
					byte c=buf[i];
					i++;
					if(c<=slashr){mode++; break;}
					if(c!=at){hdr.append((char)c);}
				}
				while(i<lim && buf[i]<=slashr){i++;} //In case of \n\r
			}
			
			if(mode==1){
				while(i<lim){
					byte c=buf[i];
					i++;
					if(c<=slashr){
						mode++;
						if(hdr.length()>0 || shared.Vector.sum(counts)>0){
							if(tsw!=null){tsw.print(toString2(hdr, counts));}else if(!SUMMARY_ONLY){System.out.print(toString2(hdr, counts));}
							hdr.setLength(0);
							for(int j=0; j<counts.length; j++){
								overall[j]+=counts[j];
								counts[j]=0;
							}
						}
						break;
					}
					counts[charToNum[c]]++;
				}
				while(i<lim && buf[i]<=slashr){i++;} //In case of \n\r
			}
			
			if(mode==2){
				while(i<lim){
					byte c=buf[i];
					i++;
					if(c<=slashr){mode++; break;}
				}
				while(i<lim && buf[i]<=slashr){i++;} //In case of \n\r
			}
			
			if(mode==3){
				while(i<lim){
					byte c=buf[i];
					i++;
					if(c<=slashr){mode=0; break;}
				}
				while(i<lim && buf[i]<=slashr){i++;} //In case of \n\r
			}
			
			if(i>=lim){
				i=0;
				lim=is.read(buf);
				limsum+=lim;
			}
		}
		
		if(hdr.length()>0 || shared.Vector.sum(counts)>0){
			if(tsw!=null){tsw.print(toString2(hdr, counts));}else if(!SUMMARY_ONLY){System.out.println(toString2(hdr, counts));}
			hdr.setLength(0);
			for(int j=0; j<counts.length; j++){
				overall[j]+=counts[j];
				counts[j]=0;
			}
		}
		
		if(tsw!=null){
			tsw.poisonAndWait();
		}
		LIMSUM=limsum;
		return overall;
	}
	
	private static String toString2(StringBuilder sb, int[] counts){
		final long sum1=(long)counts[0]+(long)counts[1]+(long)counts[2]+(long)counts[3];
		final long sum2=sum1+counts[4];
		final float inv1=1f/Tools.max(1, sum1);
		final float inv2=1f/Tools.max(1, sum2);
		if(FORMAT==1){
			return sb.append(Tools.format("\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
					sum2, counts[0]*inv1, counts[1]*inv1, counts[2]*inv1, counts[3]*inv1, counts[4]*inv2, (counts[1]+counts[2])*inv1)).toString();
		}else if(FORMAT==2){
			return sb.append(Tools.format("\t%.5f\n", (counts[1]+counts[2])*inv1)).toString();
		}else if(FORMAT==4){
			return sb.append(Tools.format("\t%d\t%.5f\n", sum2, (counts[1]+counts[2])*inv1)).toString();
		}else{
			throw new RuntimeException("Unknown format.");
		}
	}
	
	private static String toString2(StringBuilder sb, long[] counts){
		final long sum1=(long)counts[0]+(long)counts[1]+(long)counts[2]+(long)counts[3];
		final long sum2=sum1+counts[4];
		final float inv1=1f/Tools.max(1, sum1);
		final float inv2=1f/Tools.max(1, sum2);
		if(FORMAT==1){
			return sb.append(Tools.format("\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
					sum2, counts[0]*inv1, counts[1]*inv1, counts[2]*inv1, counts[3]*inv1, counts[4]*inv2, (counts[1]+counts[2])*inv1)).toString();
		}else if(FORMAT==2){
			return sb.append(Tools.format("\t%.5f\n", (counts[1]+counts[2])*inv1)).toString();
		}else if(FORMAT==4){
			return sb.append(Tools.format("\t%d\t%.5f\n", sum2, (counts[1]+counts[2])*inv1)).toString();
		}else{
			throw new RuntimeException("Unknown format.");
		}
	}
	
	private static final byte[] charToNum=makeCharToNum();
	public static int FORMAT=1;
	public static boolean SUMMARY_ONLY=false;
	private static long LIMSUM=0;

	final static byte slashr='\r', slashn='\n', carrot='>', at='@';
	
	static PrintStream outstream=System.err;
	
	/**
	 * @return
	 */
	private static byte[] makeCharToNum() {
		byte[] r=new byte[256];
		Arrays.fill(r, (byte)4);
		r['a']=r['A']=0;
		r['c']=r['C']=1;
		r['g']=r['G']=2;
		r['t']=r['T']=3;
		r['\n']=r['\r']=r['>']=r['@']=r['+']=5;
		return r;
	}
}
