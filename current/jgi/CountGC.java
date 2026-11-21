package jgi;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Counts base composition in FASTA or unwrapped FASTQ with a serial byte parser.
 * Emits per-record rows and an overall row in the selected format. A/C/G/T and GC
 * fractions use the defined-base count; the unknown-base fraction uses total bases.
 * Other ASCII sequence characters generally enter the unknown bucket; selected
 * separator symbols are ignored by the formatted base totals. These parsers are
 * not general FASTA/FASTQ validators. Record and overall buckets are signed longs.
 * The main method supports managed compressed inputs and borrowed standard input;
 * direct counting methods borrow an already-decoded input stream.
 * Global output configuration and byte accounting persist between calls; main also
 * changes shared parser/decompression settings. Concurrent use is unsafe.
 *
 * @author Brian Bushnell
 * @date Dec 13, 2012
 */
public class CountGC{
	
	/**
	 * Parses legacy arguments, counts one input and prints summary/timing diagnostics.
	 * No arguments prints help; otherwise an absent input is rejected. Formats 1, 2
	 * and 4 select the same row layout for per-record and overall output. Without a
	 * named output, SUMMARY_ONLY suppresses per-record rows. Benchmark mode consumes
	 * bytes without parsing sequences. Input/output aliases
	 * and decoder ownership follow the managed-input selection below.
	 * @param args Launcher arguments; unknown arguments are currently ignored
	 * @throws IOException If input reading or an owned input/output reports an I/O failure
	 */
	public static void main(String[] args) throws IOException{
		
		Timer t=new Timer();
		if(args.length==0){
			System.out.println("Usage: CountGC in=<infile> out=<outfile>");
			System.out.println("Alternately, 'out=stdout' will print to standard out.");
			System.out.println("Optional flag, format:");
			System.out.println("format=1\tid length A C G T N GC");
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
				final String[] split=Parse.splitOnFirst(arg, '=');
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
		Tools.testForDuplicateFiles(true, in, out);
		
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
		
		// Literal stdin/standardin is borrowed. Other inputs use an owned wrapper,
		// which keeps borrowed System.in open when ReadWrite selects a stdin alias.
		// Closing the wrapper also checks helper-tracked decoder completion; see
		// ReadWrite.getManagedInputStream for remaining backend/setup limitations.
		// The FileFormat prefix probe is completed before this parser stream opens.
		final FileFormat ff=benchmark ? null :
			FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
		if(is==null){is=ReadWrite.getManagedInputStream(in, false, true);}
		Throwable primary=null;
		try{
			if(benchmark){sum=bench2(is);}
			else{
				boolean fastq=ff.fastq();
				boolean fasta=!fastq; //Default.
				if(fastq){counts=countFastq(is, out);}
				else if(fasta){counts=countFasta(is, out);}
				else{throw new RuntimeException("Unknown or unsupported file format.");}
			}
		}catch(IOException e){
			primary=e;
			throw e;
		}catch(RuntimeException e){
			primary=e;
			throw e;
		}catch(Error e){
			primary=e;
			throw e;
		}finally{
			closeOwnedInput(is, primary);
		}
		

		t.stop();
		
		if(benchmark){
			outstream.println("Time: \t"+t);
			long bytes=new File(in).length();
			if(bytes<1){bytes=sum;}
			double mbps1=bytes*1000d/t.elapsed;
			double mbps2=sum*1000d/t.elapsed;
			outstream.println(Tools.format("Raw Speed:         \t%.2f MBytes/s", mbps1));
			outstream.println(Tools.format("Uncompressed Speed:\t%.2f MBytes/s", mbps2));
		}else{
			outstream.println(toString2(new StringBuilder("Overall"), counts));
			outstream.println("Time: \t"+t);
			long bytes=new File(in).length();
			if(bytes<1){bytes=LIMSUM;}
			double mbps=bytes*1000d/t.elapsed;
			double mbpps=simd.Vector.sum(counts)*1000d/t.elapsed;
			outstream.println(Tools.format("Speed:\t%.2f MBytes/s", mbps));
			outstream.println(Tools.format("      \t%.2f MBases/s", mbpps));
		}
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}

	/**
	 * Closes CountGC.main's owned non-System.in stream, suppressing close
	 * failures onto any existing parser/read failure.
	 */
	private static void closeOwnedInput(final InputStream is, final Throwable primary)
			throws IOException{
		if(is==null || is==System.in){return;}
		try{
			is.close();
		}catch(IOException e){
			if(primary==null){throw e;}
			else if(primary!=e){primary.addSuppressed(e);}
		}catch(RuntimeException e){
			if(primary==null){throw e;}
			else if(primary!=e){primary.addSuppressed(e);}
		}catch(Error e){
			if(primary==null){throw e;}
			else if(primary!=e){primary.addSuppressed(e);}
		}
	}
	
	/**
	 * Consumes a borrowed stream in 32-KiB chunks until read returns zero or EOF.
	 * Does not parse records, emit rows, close the stream or update LIMSUM.
	 * @param is Stream already positioned at the first byte to consume
	 * @return Number of bytes returned by positive reads
	 * @throws IOException If reading fails
	 */
	public static long bench2(InputStream is) throws IOException{
		final byte[] buf=new byte[32768];
		long sum=0;
		for(long len=is.read(buf); len>0; len=is.read(buf)){sum+=len;}
		return sum;
	}
	
	/**
	 * Counts base composition in FASTA format sequences.
	 * Processes sequences one at a time with long record and overall buckets.
	 * A greater-than sign starts a header; ASCII values at most CR terminate headers.
	 * Non-ASCII sequence bytes fail, while header bytes are decoded as UTF-8 for output.
	 * Outputs per-sequence statistics and returns overall counts.
	 * Leaves the caller-owned input stream open. Flushes and joins the owned
	 * output writer on normal completion and parser/read failure, preserving
	 * any primary failure.
	 *
	 * @param is Input stream containing FASTA sequences
	 * @param out Output file path, or null for stdout/summary only
	 * @return Six overall buckets [A, C, G, T, unknown ASCII, ignored separator symbols]
	 * @throws IOException If stream reading fails, or an output write/closure
	 * failure is reported
	 */
	public static long[] countFasta(InputStream is, String out) throws IOException{
		
		long limsum=0;
		final byte[] buf=new byte[32768];
		final TextStreamWriter tsw=(out==null ? null : new TextStreamWriter(out, true, false, false));
		if(tsw!=null){tsw.start();}
		Throwable primary=null;
		try{
			final long[] counts=KillSwitch.allocLong1D(6);
			final long[] overall=new long[6];
			final ByteBuilder hdr=new ByteBuilder();
			boolean hdmode=false;

			int i=0;
			int lim=is.read(buf);
			if(lim>0){limsum+=lim;}
			
			while(lim>0){
				if(hdmode){
					while(i<lim){
						byte c=buf[i];
						i++;
						if(isLineBreakOrControl(c)){hdmode=false; break;}
						hdr.append(c);
					}
				}

				if(!hdmode){
					while(i<lim){
						byte c=buf[i];
						i++;

						if(c==carrot){
							hdmode=true;
							if(hdr.length()>0 || simd.Vector.sum(counts)>0){
								if(tsw!=null){tsw.print(toString2(hdr, counts));}else if(!SUMMARY_ONLY){System.out.print(toString2(hdr, counts));}
								hdr.clear();
								for(int j=0; j<counts.length; j++){
									overall[j]+=counts[j];
									counts[j]=0;
								}
							}
							break;
						}
						counts[baseToNumber(c, "FASTA")]++;

					}
				}
				if(i>=lim){
					i=0;
					lim=is.read(buf);
					if(lim>0){limsum+=lim;}
				}
			}

			if(hdr.length()>0 || simd.Vector.sum(counts)>0){
				if(tsw!=null){tsw.print(toString2(hdr, counts));}else if(!SUMMARY_ONLY){System.out.println(toString2(hdr, counts));}
				hdr.clear();
				for(int j=0; j<counts.length; j++){
					overall[j]+=counts[j];
					counts[j]=0;
				}
			}

			LIMSUM=limsum;
			return overall;
		}catch(IOException e){
			primary=e;
			throw e;
		}catch(RuntimeException e){
			primary=e;
			throw e;
		}catch(Error e){
			primary=e;
			throw e;
		}finally{
			finishWriter(tsw, primary);
		}
	}
	
	/**
	 * Counts base composition in FASTQ format sequences.
	 * Expects four-line records and counts only sequence lines; this is not a strict
	 * framing validator and does not check plus markers or sequence/quality lengths.
	 * Runs of ASCII values at most CR separate lines, including across buffer refills.
	 * Removes only the first header at-sign and decodes retained header bytes as UTF-8.
	 * Record and accumulated overall buckets are longs.
	 * Outputs per-sequence statistics and returns overall counts.
	 * Leaves the caller-owned input stream open. Flushes and joins the owned
	 * output writer on normal completion and parser/read failure, preserving
	 * any primary failure.
	 *
	 * @param is Input stream containing FASTQ sequences
	 * @param out Output file path, or null for stdout/summary only
	 * @return Six overall buckets [A, C, G, T, unknown ASCII, ignored separator symbols]
	 * @throws IOException If stream reading fails, or an output write/closure
	 * failure is reported
	 */
	public static long[] countFastq(InputStream is, String out) throws IOException{
//		assert(false) : "Fastq mode - TODO"; //TODO
		long limsum=0;
		final byte[] buf=new byte[32768];
		final TextStreamWriter tsw=(out==null ? null : new TextStreamWriter(out, true, false, false));
		if(tsw!=null){tsw.start();}
		Throwable primary=null;
		try{
			final long[] counts=KillSwitch.allocLong1D(6);
			final long[] overall=new long[6];
			final ByteBuilder hdr=new ByteBuilder();

			int mode=3;
			boolean firstHeaderByte=true;

			int i=0;
			int lim=is.read(buf);
			if(lim>0){limsum+=lim;}

			while(mode==3 && lim>0){
				while(i<lim && buf[i]!=at){i++;}
				if(i>=lim){
					lim=is.read(buf);
					if(lim>0){limsum+=lim;}
				}else{
					assert(buf[i]==at);
					mode=0;
				}
			}
			
			// Carries a line-break/control run across buffer refills.
			boolean skipBreaks=false;
			while(lim>0){
				if(skipBreaks){
					while(i<lim && isLineBreakOrControl(buf[i])){i++;}
					if(i>=lim){
						i=0;
						lim=is.read(buf);
						if(lim>0){limsum+=lim;}
						continue;
					}
					skipBreaks=false;
				}

				if(mode==0){
					boolean endedLine=false;
					while(i<lim){
						byte c=buf[i];
						i++;
						if(isLineBreakOrControl(c)){
							mode++;
							firstHeaderByte=true;
							endedLine=true;
							break;
						}
						if(!firstHeaderByte || c!=at){hdr.append(c);}
						firstHeaderByte=false;
					}
					while(i<lim && isLineBreakOrControl(buf[i])){i++;} //In case of \n\r
					if(endedLine && i>=lim){skipBreaks=true;}
				}

				if(mode==1){
					boolean endedLine=false;
					while(i<lim){
						byte c=buf[i];
						i++;
						if(isLineBreakOrControl(c)){
							mode++;
							endedLine=true;
							if(hdr.length()>0 || simd.Vector.sum(counts)>0){
								if(tsw!=null){tsw.print(toString2(hdr, counts));}else if(!SUMMARY_ONLY){System.out.print(toString2(hdr, counts));}
								hdr.clear();
								for(int j=0; j<counts.length; j++){
									overall[j]+=counts[j];
									counts[j]=0;
								}
							}
							break;
						}
						counts[baseToNumber(c, "FASTQ")]++;
					}
					while(i<lim && isLineBreakOrControl(buf[i])){i++;} //In case of \n\r
					if(endedLine && i>=lim){skipBreaks=true;}
				}

				if(mode==2){
					boolean endedLine=false;
					while(i<lim){
						byte c=buf[i];
						i++;
						if(isLineBreakOrControl(c)){
							mode++;
							endedLine=true;
							break;
						}
					}
					while(i<lim && isLineBreakOrControl(buf[i])){i++;} //In case of \n\r
					if(endedLine && i>=lim){skipBreaks=true;}
				}

				if(mode==3){
					boolean endedLine=false;
					while(i<lim){
						byte c=buf[i];
						i++;
						if(isLineBreakOrControl(c)){
							mode=0;
							endedLine=true;
							break;
						}
					}
					while(i<lim && isLineBreakOrControl(buf[i])){i++;} //In case of \n\r
					if(endedLine && i>=lim){skipBreaks=true;}
				}

				if(i>=lim){
					i=0;
					lim=is.read(buf);
					if(lim>0){limsum+=lim;}
				}
			}
			
			if(hdr.length()>0 || simd.Vector.sum(counts)>0){
				if(tsw!=null){tsw.print(toString2(hdr, counts));}else if(!SUMMARY_ONLY){System.out.println(toString2(hdr, counts));}
				hdr.clear();
				for(int j=0; j<counts.length; j++){
					overall[j]+=counts[j];
					counts[j]=0;
				}
			}
			
			LIMSUM=limsum;
			return overall;
		}catch(IOException e){
			primary=e;
			throw e;
		}catch(RuntimeException e){
			primary=e;
			throw e;
		}catch(Error e){
			primary=e;
			throw e;
		}finally{
			finishWriter(tsw, primary);
		}
	}

	/** Throws output IOException directly or suppresses it onto existing parser/read failures. */
	private static void finishWriter(final TextStreamWriter tsw, final Throwable primary)
			throws IOException{
		if(tsw==null){return;}
		try{
			final boolean error=tsw.poisonAndWait();
			if(error){
				final IOException ioe=new IOException("CountGC output error for "+tsw.fname);
				if(primary!=null){primary.addSuppressed(ioe);}
				else{throw ioe;}
			}
		}catch(RuntimeException e){
			if(primary!=null){primary.addSuppressed(e);}
			else{throw e;}
		}catch(Error e){
			if(primary!=null){primary.addSuppressed(e);}
			else{throw e;}
		}
	}
	
	/** Decodes the live UTF-8 header bytes, then appends a formatted count row. */
	private static String toString2(ByteBuilder hdr, long[] counts){
		final String name=new String(hdr.array, 0, hdr.length(), StandardCharsets.UTF_8);
		return toString2(new StringBuilder(name), counts);
	}

	/**
	 * Appends one FORMAT-selected row and its newline to the caller's header builder.
	 * Excludes bucket 5; fractions use float reciprocals with denominators clamped to one.
	 * Widens record-bucket sums to long but does not repair overflow within a bucket.
	 */
	private static String toString2(StringBuilder sb, int[] counts){
		final long sum1=(long)counts[0]+(long)counts[1]+(long)counts[2]+(long)counts[3];
		final long sum2=sum1+counts[4];
		final long gc=(long)counts[1]+(long)counts[2];
		final float inv1=1f/Tools.max(1, sum1);
		final float inv2=1f/Tools.max(1, sum2);
		if(FORMAT==1){
			return sb.append(Tools.format("\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
					sum2, counts[0]*inv1, counts[1]*inv1, counts[2]*inv1, counts[3]*inv1, counts[4]*inv2, gc*inv1)).toString();
		}else if(FORMAT==2){
			return sb.append(Tools.format("\t%.5f\n", gc*inv1)).toString();
		}else if(FORMAT==4){
			return sb.append(Tools.format("\t%d\t%.5f\n", sum2, gc*inv1)).toString();
		}else{
			throw new RuntimeException("Unknown format.");
		}
	}
	
	/** Appends a record or overall row; long sums are unchecked and bucket 5 is excluded. */
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
	
	/** ASCII-to-bucket lookup; negative signed sequence bytes are rejected by baseToNumber. */
	private static final byte[] charToNum=makeCharToNum();
	/** Output row layout: 1 for composition, 2 for GC alone, 4 for length and GC. */
	public static int FORMAT=1;
	/** Suppresses per-record stdout rows when no named output writer exists. */
	public static boolean SUMMARY_ONLY=false;
	/** Sum of positive input-read lengths, assigned after parsing and before output-writer finalization. */
	private static long LIMSUM=0;

	// ASCII line breaks and FASTA/FASTQ header prefix bytes.
	final static byte slashr='\r', slashn='\n', carrot='>', at='@';
	
	/** Main's summary/timing stream, replaceable through PreParser; not the per-record stdout path. */
	static PrintStream outstream=System.err;
	
	/** Returns a base-table value or throws a clear error for non-ASCII sequence bytes. */
	private static int baseToNumber(final byte c, final String format){
		try{
			return charToNum[c];
		}catch(ArrayIndexOutOfBoundsException e){
			// Fixed 256-slot table makes signed high-bit bytes the only misses.
			throw new RuntimeException("Invalid "+format+" sequence byte "+toHex(c));
		}
	}

	/** Tests whether a byte is an ASCII line break or lower control character. */
	private static boolean isLineBreakOrControl(final byte c){
		return (c&0xFF)<=slashr;
	}

	/** Formats a byte as a two-digit hexadecimal literal. */
	private static String toHex(final byte c){
		final String s=Integer.toHexString(c&0xFF).toUpperCase();
		return s.length()==1 ? "0x0"+s : "0x"+s;
	}

	/**
	 * Creates lookup table for converting ASCII characters to base indices.
	 * Maps bytes 0..127 as A/a=0, C/c=1, G/g=2, T/t=3 and other ASCII=4,
	 * except CR, LF, greater-than, at-sign and plus, which map to ignored bucket 5.
	 * Slots 128..255 contain -1; signed-byte indexing rejects such sequence bytes earlier.
	 * @return Byte array mapping ASCII values to base type indices
	 */
	private static byte[] makeCharToNum(){
		byte[] r=new byte[256];
		Arrays.fill(r, (byte)-1);
		Arrays.fill(r, 0, 128, (byte)4);
		r['a']=r['A']=0;
		r['c']=r['C']=1;
		r['g']=r['G']=2;
		r['t']=r['T']=3;
		r['\n']=r['\r']=r['>']=r['@']=r['+']=5;
		return r;
	}
}
