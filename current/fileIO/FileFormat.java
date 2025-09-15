package fileIO;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jgi.TestFormat;
import shared.Parse;
import shared.PreParser;
import shared.Tools;
import stream.Read;

/**
 * This class contains metadata about a file
 * @author Brian Bushnell
 * @date Dec 19, 2012
 *
 */
public final class FileFormat {
	
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null /*new Object() { }.getClass().getEnclosingClass()*/, false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		stream.FASTQ.warnQualityChange=false;
		PRINT_WARNING=false;
		boolean full=false;
		ArrayList<String> files=new ArrayList<String>();
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("full")){
				full=Parse.parseBoolean(b);
			}else if(b!=null){
//				assert(a.startsWith("in")) : "Unknown parameter "+arg;
				if(a.startsWith("in")){files.add(b);}
			}else{
				files.add(arg);
			}
		}
		
		if(full){
			TestFormat.main(args);
		}else{
			for(String fname : files){
				test(fname, true);
			}
		}
		
	}
	
	private static void test(String fname, boolean forceFileRead){
		FileFormat ffName=testInput(fname, FASTQ, null, false, false, false);
		FileFormat ffContent=testInput(fname, ffName.format(), null, false, true, true);
		FileFormat ff=ffContent;
//		assert(false) : ffName+"\n"+ffContent;
		if(ff==null){
			System.out.println("null");
		}else{
			int q=33;
			int len=-1;
			boolean i=false;
			if(ff.fastq()){
				byte qold=stream.FASTQ.ASCII_OFFSET;
				stream.FASTQ.ASCII_OFFSET=33;
				int[] qi=testInterleavedAndQuality(fname, false);
				q=qi[0];
				i=(qi[1]==INTERLEAVED);
				len=qi[2];
				stream.FASTQ.ASCII_OFFSET=qold;
				//TODO: Check barcode status here
			}else if(ff.fasta()){
				i=stream.FASTQ.testInterleavedFasta(fname, false);
			}
			if(ff.isSequence()){
				String qs=(q==33 ? "sanger" : q==64 ? "illumina" : ""+q);
				System.out.print(qs+"\t"+FORMAT_ARRAY[ff.format()]+"\t"+COMPRESSION_ARRAY[ff.compression()]);
				System.out.print("\t"+(i ? "interleaved" : "single-ended"));
				if(len>0){System.out.print("\t"+len+"bp");}
			}else{
				System.out.print(FORMAT_ARRAY[ff.format()]+"\t"+COMPRESSION_ARRAY[ff.compression()]);
			}
			if(ffName.format()!=ff.format()){System.out.print("\t"+FORMAT_ARRAY[ffName.format()]+"\t(File extension differs from contents)");}
			System.out.println();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static FileFormat testInput(String fname, String overrideExtension, boolean allowSubprocess){
		if(verbose){System.err.println("testInputA("+fname+", "+overrideExtension+", "+allowSubprocess+")");}
		return testInput(fname, FASTQ, overrideExtension, allowSubprocess, true);
	}
	
	public static FileFormat[] testInputList(List<String> fname, int defaultFormat, String overrideExtension, boolean allowSubprocess, boolean allowFileRead){
		if(verbose){System.err.println("testInputList("+fname+", "+defaultFormat+", "+overrideExtension+", "+allowSubprocess+", "+allowFileRead+")");}
		FileFormat[] ffa=new FileFormat[fname.size()];
		for(int i=0; i<fname.size(); i++){
			ffa[i]=testInput(fname.get(i), defaultFormat, overrideExtension, allowSubprocess, allowFileRead, false);
		}
		return ffa;
	}
	
	public static FileFormat[] testInput(String fnames[], int defaultFormat, String overrideExtension, boolean allowSubprocess, boolean allowFileRead){
		FileFormat[] array=new FileFormat[fnames.length];
		for(int i=0; i<fnames.length; i++){
			array[i]=testInput(fnames[i], defaultFormat, overrideExtension, allowSubprocess, allowFileRead);
		}
		return array;
	}
	
	public static FileFormat testInput(String fname, int defaultFormat, String overrideExtension, boolean allowSubprocess, boolean allowFileRead){
		if(verbose){System.err.println("testInputB("+fname+", "+defaultFormat+", "+overrideExtension+", "+allowSubprocess+", "+allowFileRead+")");}
		return testInput(fname, defaultFormat, overrideExtension, allowSubprocess, allowFileRead, allowFileRead && !isStdin(fname));
	}
	
	public static FileFormat testInput(String fname, int defaultFormat, String overrideExtension, boolean allowSubprocess, boolean allowFileRead, boolean forceFileRead){
		if(verbose){System.err.println("testInputC("+fname+", "+defaultFormat+", "+overrideExtension+", "+allowSubprocess+", "+allowFileRead+", "+forceFileRead+")");}
		if(fname==null){return null;}
		int overrideFormat=0;
		int overrideCompression=0;
		if(overrideExtension!=null && overrideExtension.length()>0){
			int[] a=testFormat(overrideExtension, false, false);
			if(a!=null){
				overrideFormat=a[0];
				if(a[1]!=RAW){overrideCompression=a[1];}
			}
		}
		return testInput(fname, defaultFormat, overrideFormat, overrideCompression, allowSubprocess, allowFileRead, forceFileRead);
	}
	
	/**
	 * Create an input FileFormat object for this filename.
	 * @param fname Filename (path).
	 * @param defaultFormat Use this format if the name is unclear and the format is not autodetected.
	 * @param overrideFormat If specified, ignore the file extension and autodetection and input using this format.
	 * @param overrideCompression If specified, ignore the file extension and input using this compression protocol.
	 * @param allowSubprocess Permission to spawn a subprocess like bgzip.
	 * @param allowFileRead Permission to read the file while constructing this FileFormat, for the purpose of format detection.
	 * @param forceFileRead Force reading the file while constructing this FileFormat, for the purpose of format detection.
	 * @return A FileFormat, or null if the filename is null.
	 */
	public static FileFormat testInput(String fname, int defaultFormat, int overrideFormat,
			int overrideCompression, boolean allowSubprocess, boolean allowFileRead, boolean forceFileRead){
		if(verbose){System.err.println("testInputD("+fname+", "+defaultFormat+", "+overrideFormat+", "+overrideCompression+", "+allowSubprocess+", "+allowFileRead+", "+forceFileRead+")");}
		if(fname==null){return null;}
		return new FileFormat(fname, READ, defaultFormat, overrideFormat, overrideCompression, allowSubprocess, allowFileRead, forceFileRead, false, false, false, true);
	}
	
	/**
	 * Create an output FileFormat object for this filename.
	 * @param fname Filename (path).
	 * @param defaultFormat Use this format if the name is unclear.
	 * @param overrideExtension If specified, ignore the file extension and output in this format.
	 * @param allowSubprocess Permission to spawn a subprocess like bgzip.
	 * @param overwrite Permission to overwrite existing files.
	 * @param append Permission to append to existing files.
	 * @param ordered True if the input order should be maintained (for multithreaded read processing).
	 * @return A FileFormat, or null if the filename is null.
	 */
	public static FileFormat testOutput(String fname, int defaultFormat, String overrideExtension, boolean allowSubprocess, boolean overwrite, boolean append, boolean ordered){
		if(fname==null){return null;}
		int overrideFormat=0;
		int overrideCompression=0;
		if(overrideExtension!=null && overrideExtension.length()>0){
			int[] a=testFormat(overrideExtension, false, false);
			if(a!=null){
				overrideFormat=a[0];
				if(a[1]!=RAW){overrideCompression=a[1];}
			}
		}
		return testOutput(fname, defaultFormat, overrideFormat, overrideCompression, allowSubprocess, overwrite, append, ordered);
	}
	
	public static FileFormat testOutput(String fname, int defaultFormat, int overrideFormat, int overrideCompression, boolean allowSubprocess, boolean overwrite, boolean append, boolean ordered){
		if(fname==null){return null;}
		return new FileFormat(fname, WRITE, defaultFormat, overrideFormat, overrideCompression, allowSubprocess, false, false, overwrite, append, ordered, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Constructor         ----------------*/
	/*--------------------------------------------------------------*/
	
	private FileFormat(String fname, int mode_, int defaultFormat, int overrideFormat, int overrideCompression, boolean allowSubprocess_,
			boolean allowFileRead, boolean forceFileRead, boolean overwrite_, boolean append_, boolean ordered_, boolean input_){
//			, boolean interleaved_, long maxReads_){
		
		if(verbose){
//			new Exception().printStackTrace(System.err);
			System.err.println("FileFormat(fname="+fname+", mode="+mode_+", dFormat="+defaultFormat+", oFormat="+overrideFormat+", oCompression="+overrideCompression+
					", allowSub="+allowSubprocess_+", allowRead="+allowFileRead+", forceFileRead="+forceFileRead+
					", ow="+overwrite_+", append="+append_+", ordered="+ordered_+")");
		}
		assert(!forceFileRead || allowFileRead);
		
//		assert(!overwrite_ || !append_) : "Both overwrite and append may not be set to true.";
		if(overwrite_ && append_){overwrite_=false;}
		
		assert(fname!=null);
		fname=fname.trim().replace('\\', '/');
		assert(fname.trim().length()>0) : fname;
		
		if(defaultFormat<1 && !forceFileRead){defaultFormat=FQ;}
		allowFileRead&=(mode_==READ);
		int[] a=testFormat(fname, allowFileRead, forceFileRead);
		
		if(verbose){System.err.println(Arrays.toString(a));}
		
		if(a[0]==UNKNOWN && overrideFormat<1){
			a[0]=defaultFormat;
			if(defaultFormat!=TEXT && PRINT_WARNING && a[2]!=DEVNULL){
				System.err.println("Unspecified format for "+(mode_==READ ? "input" : "output")+
						" "+(fname==null ? "stream" : fname)+"; defaulting to "+FORMAT_ARRAY[a[0]]+".");
			}
		}
		if(verbose){System.err.println(Arrays.toString(a));}
		
		if(overrideFormat>0){a[0]=overrideFormat;}
		if(overrideCompression>0){a[1]=overrideCompression;}
		
		if(verbose){System.err.println(Arrays.toString(a));}

		
//		{format, compression, type, interleaved, quality, length}
		name=fname;
		simpleName=new File(name).getName();
		format=a[0];
		compression=a[1];
		type=a[2];
		interleaving=a[3];
		asciiOffset=a[4];
		length=a[5];
		numBarcodes=a[6];
		barcodeDelimiter=a[7];
		barcodeLength1=a[8];
		barcodeLength2=a[9];
		mode=mode_;
		input=input_;
//		assert(!allowFileRead || barcodeDelimiter>=0) : barcodeDelimiter;
		
		overwrite=overwrite_;
		append=append_;
		allowSubprocess=allowSubprocess_;
		ordered=ordered_;
		amino=isAminoExt(rawExtension());
		
//		interleaved=interleaved_;
//		maxReads=write() ? -1 : maxReads_;

		assert(forceFileRead || !unknownFormat()) : "Unknown file format for "+fname+"\n"+
			mode_+", "+defaultFormat+", "+overrideFormat+", "+overrideCompression+", "+allowSubprocess_+", "+allowFileRead+", "+overwrite_;
		assert(!unknownCompression()) : "Unknown compression for "+fname+"\n"+
			mode_+", "+defaultFormat+", "+overrideFormat+", "+overrideCompression+", "+allowSubprocess_+", "+allowFileRead+", "+overwrite_;
		assert(!unknownType()) : "Unknown stream type for "+fname+"\n"+
			mode_+", "+defaultFormat+", "+overrideFormat+", "+overrideCompression+", "+allowSubprocess_+", "+allowFileRead+", "+overwrite_;
		assert(!unknownMode()) : "Unknown I/O mode for "+fname+"\n"+
			mode_+", "+defaultFormat+", "+overrideFormat+", "+overrideCompression+", "+allowSubprocess_+", "+allowFileRead+", "+overwrite_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append(name).append(',');
		sb.append(format+"("+FORMAT_ARRAY[format]+")").append(',');
		sb.append(compression+"("+COMPRESSION_ARRAY[compression]+")").append(',');
		sb.append(type+"("+TYPE_ARRAY[type]+")").append(',');
		sb.append(interleaving+"("+INTERLEAVING_ARRAY[interleaving]+")").append(',');
//		sb.append("ascii"+asciiOffset).append(',');
		sb.append(mode+"("+MODE_ARRAY[mode]+")").append(',');
		sb.append("ow="+(overwrite ? "t" : "f")).append(',');
		sb.append("app="+(append ? "t" : "f")).append(',');
		sb.append("sub="+(allowSubprocess ? "t" : "f")).append(',');
		sb.append("ordered="+(ordered ? "t" : "f"));
		return sb.toString();
	}
	
	public static String toString(int[] vector){
		/* {format, compression, type, interleaved, 
		 * quality, length, numBarcodes, barcodeDelimiter} */
		
		int format=vector[0], compression=vector[1], type=vector[2], interleaving=vector[3];
		int quality=vector[4], length=vector[5], codes=vector[6], delimiter=vector[7];
		StringBuilder sb=new StringBuilder();
		sb.append(format+"("+FORMAT_ARRAY[format]+")").append(',');
		sb.append(compression+"("+COMPRESSION_ARRAY[compression]+")").append(',');
		sb.append(type+"("+TYPE_ARRAY[type]+")").append(',');
		sb.append(interleaving+"("+INTERLEAVING_ARRAY[interleaving]+")").append(',');

		sb.append("ASCII-"+quality).append(',');
		sb.append("length "+length).append(',');
		sb.append("barcodes "+codes).append(',');
		sb.append("delimiter "+(delimiter>0 ? Character.toString(delimiter) : Integer.toString(delimiter)));
		
		
		return sb.toString();
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns an int array: {format, compression, type, interleaved, 
	 * quality, length, numBarcodes, barcodeDelimiter, bcLength1, bcLength2} */
	public static final int[] testFormat(String fname, boolean allowFileRead, boolean forceFileRead){
		if(verbose){System.err.println("testFormat("+fname+", "+allowFileRead+", "+forceFileRead+")");}
		final int[] r=new int[] {UNKNOWN, RAW, FILE, UNKNOWN, -1, -1, -1, -1, -1, -1};
		if(fname==null || fname.length()<1){
			r[2]=STDIO;
			return r;
		}
		String slc=fname.trim().toLowerCase();
		if(slc.indexOf('/')<0){slc=slc.substring(slc.lastIndexOf('/')+1);}
		String comp=ReadWrite.compressionType(slc);
		String ext=ReadWrite.rawExtension(slc);
		
		if(ext==null){}
		else if(ext.equals("fq") || ext.equals("fastq") || (comp!=null && comp.equals("fqz"))){r[0]=FASTQ;}
		else if(isFastaExt(ext)){r[0]=FASTA;}
		else if(ext.equals("bread")){r[0]=BREAD;}
		else if(ext.equals("sam")){r[0]=SAM;}
		else if(ext.equals("csfasta")){r[0]=CSFASTA;}
		else if(ext.equals("qual")){r[0]=QUAL;}
		else if(ext.equals("bam")){r[0]=BAM;}
		else if(ext.equals("sites") || ext.equals("sitesonly")){r[0]=SITES;}
		else if(ext.equals("info") || ext.equals("attachment")){r[0]=ATTACHMENT;}
		else if(ext.equals("scarf")){r[0]=SCARF;}
		else if(ext.equals("phylip")){r[0]=PHYLIP;}
		else if(ext.equals("header") || ext.equals("headers")){r[0]=HEADER;}
		else if(ext.equals("int1d")){r[0]=INT1D;}
		else if(ext.equals("long1d")){r[0]=LONG1D;}
		else if(ext.equals("bitset")){r[0]=BITSET;}
		else if(ext.equals("sketch")){r[0]=SKETCH;}
		else if(ext.equals("oneline") || ext.equals("flat")){r[0]=ONELINE;}
		else if(ext.equals("fastr") || ext.equals("fr")){r[0]=FASTR;}
		else if(ext.equals("vcf")){r[0]=VCF;}
		else if(ext.equals("var")){r[0]=VAR;}
		else if(ext.equals("gff") || ext.equals("gff3")){r[0]=GFF;}
		else if(ext.equals("bed")){r[0]=BED;}
		else if(ext.equals("pgm") || ext.equals("pkm")){r[0]=PGM;}
		else if(ext.equals("embl")){r[0]=EMBL;}
		else if(ext.equals("gbk")){r[0]=GBK;}
		else if(ext.equals("gbff")){r[0]=GBFF;}
		else if(ext.equals("alm")){r[0]=ALM;}
		else if(ext.equals("bbnet")){r[0]=BBNET;}
		else if(ext.equals("bbvec") || ext.equals("vec")){r[0]=BBVEC;}
		else if(ext.equals("clade") || ext.equals("spectra")){r[0]=CLADE;}
		
		if(comp!=null){
			r[1]=Tools.find(comp, COMPRESSION_ARRAY);
			assert(r[1]>0) : "Unhandled compression type: "+comp;
		}
		
//		if(r[1]==GZIP && allowFileRead){
//			//Check magic number, perhaps
//		}
		
		if(slc.length()>2 && slc.charAt(0)=='s' && slc.charAt(1)=='t'){
			if(slc.equals("stdin") || slc.startsWith("stdin.") || slc.equals("standardin")){r[2]=STDIO;}
			else if(slc.equals("stdout") || slc.startsWith("stdout.") || slc.equals("standardout")){r[2]=STDIO;}
		}else if("/dev/null".equalsIgnoreCase(slc) || "devnull".equalsIgnoreCase(slc)){
			r[2]=DEVNULL;
		}

		if(verbose){System.err.println("Before reading: \t"+r[0]+", "+toString(r)+
				", "+forceFileRead+", "+(r[0]!=BAM));
			System.err.println((r[0]==UNKNOWN)+","+(r[0]!=BAM && forceFileRead)+","+ 
				(r[0]==FASTQ || r[0]==FASTA)+","+(r[3]==UNKNOWN)+","+allowFileRead+","
					+(!stream.FASTQ.FORCE_INTERLEAVED)
				+","+stream.FASTQ.TEST_INTERLEAVED);//false,false,true,true,true,true,false
		}
		if(r[0]==UNKNOWN || (r[0]!=BAM && forceFileRead) || 
				((r[0]==FASTQ || r[0]==FASTA) && r[3]==UNKNOWN && allowFileRead 
				&& !stream.FASTQ.FORCE_INTERLEAVED && stream.FASTQ.TEST_INTERLEAVED)){
			
			File f=(allowFileRead && r[2]==FILE ? new File(fname) : null);
			if(f!=null && f.exists() && !f.isDirectory()){
//				//a: {quality, interleaved, length, format}
//				//r: {format, compression, type, interleaved, quality, length}
				try {
					int[] a=testInterleavedAndQuality(fname, false);
					if(a!=null){
						final int aq=a[0], ai=a[1], al=a[2], af=a[3], abc=a[4], abd=a[5];
						if(aq>-1){r[4]=aq;}
						if(ai!=UNKNOWN){r[3]=ai;}
						if(r[0]==UNKNOWN || (af!=TEXT && af!=BREAD)){r[0]=af;}
						if(al>1 && r[5]==-1){r[5]=al;}
						r[6]=abc;
						r[7]=abd;
						r[8]=a[6];
						r[9]=a[7];
					}
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				if(verbose){System.err.println("After reading:   \t"+r[0]+", "+toString(r)+
						", "+forceFileRead+", "+(r[0]!=BAM));}
			}else if(r[0]==UNKNOWN){
				if(fname.equals("sequential")){r[0]=SEQUENTIAL;}
				else if(fname.equals("random")){r[0]=RANDOM;}
				else if(fname.equals("sitesonly")){r[0]=SITES;}
			}
		}
//		new Exception().printStackTrace();
		if(r[0]==FASTQ || r[0]==FASTA){
			if(r[3]==UNKNOWN) {
				if(stream.FASTQ.FORCE_INTERLEAVED){r[3]=2;}
				else{r[3]=1;}
			}else if(r[0]==FASTQ && stream.FASTQ.FORCE_INTERLEAVED) {
				r[3]=2;
			}else if(!stream.FASTQ.FORCE_INTERLEAVED && ! stream.FASTQ.TEST_INTERLEAVED) {
				r[3]=1;
			}
		}
//		assert(false) : Arrays.toString(r);
		
		if(r[2]==STDIO && allowFileRead){
			File f=new File(fname);
			if(f.exists() && !f.isDirectory()){r[2]=FILE;}
		}
		if(verbose){System.err.println("testFormat return:\t"+r[0]+", "+toString(r)+", "+forceFileRead+", "+(r[0]!=BAM)+", "+r[4]);}
		return r;
	}
	
	public static boolean hasGffExtension(String fname){
		int[] r=testFormat(fname, false, false);
		return r[0]==GFF;
	}
	
	public static boolean hasVcfExtension(String fname){
		int[] r=testFormat(fname, false, false);
		return r[0]==VCF;
	}
	
	public static boolean hasFastaExtension(String fname){
		int[] r=testFormat(fname, false, false);
		return r[0]==FA;
	}
	
	public static boolean hasFastqExtension(String fname){
		int[] r=testFormat(fname, false, false);
		return r[0]==FQ;
	}
	
	public static boolean hasFastqOrFastqExtension(String fname){
		int[] r=testFormat(fname, false, false);
		return r[0]==FQ || r[0]==FA;
	}
	
	public static boolean hasTextExtension(String fname){
		String ext=ReadWrite.rawExtension(fname);
		return ("txt".equals(ext) || "tsv".equals(ext) || "csv".equals(ext));
	}
	
	public static boolean hasSamOrBamExtension(String fname){
		int[] r=testFormat(fname, false, false);
		return r[0]==SAM || r[0]==BAM;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            ???????           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * @param fname File to read
	 * @return {quality, interleaved, length, format, numBarcodes, barcodeDelimiter,
	 * bcLen1, bcLen2}
	 */
	public static int[] testInterleavedAndQuality(String fname, boolean forceFastq){
		final ArrayList<String> oct=getFirstOctet(fname);
		return testInterleavedAndQuality(oct, fname, forceFastq);
	}
	
	public static ArrayList<String> getFirstOctet(String fname){
		if(fname==null){return null;}
		if(fname.equalsIgnoreCase("stdin") || fname.toLowerCase().startsWith("stdin.")){return null;}
		
		ArrayList<String> oct=new ArrayList<String>(8);
		
		{
			InputStream is=ReadWrite.getInputStream(fname, false, fname.toLowerCase().endsWith(".bz2"));
			BufferedReader br=new BufferedReader(new InputStreamReader(is));
			try {
				int cntr=0;
				for(String s=br.readLine(); s!=null && cntr<8; s=br.readLine()){
					oct.add(s);
					cntr++;
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			ReadWrite.finishReading(is, fname, true, br);
		}
		return oct;
	}
	
	/**
	 * @param oct First 8 lines of file
	 * @param fname File to read
	 * @return {quality, interleaved, length, format, numBarcodes, barcodeDelimiter}
	 */
	public static int[] testInterleavedAndQuality(final ArrayList<String> oct, String fname, boolean forceFastq){
		int len=-1, format=UNKNOWN;
		byte q=-1, i=UNKNOWN;
		int numBarcodes=0, barcodeDelimiter=0, bcLen1=-1, bcLen2=-1;
		if(oct==null || oct.size()<1){
			return new int[] {q, i, len, format, numBarcodes, barcodeDelimiter, bcLen1, bcLen2};
		}
		{
			String s1=oct.size()>0 ? oct.get(0) : "";
			String s2=oct.size()>1 ? oct.get(1) : "";
			String s3=oct.size()>2 ? oct.get(2) : "";
			int b1=(s1.length()>0 ? s1.charAt(0) : -1);
			int b2=(s2.length()>0 ? s2.charAt(0) : -1);
			int b3=(s3.length()>0 ? s3.charAt(0) : -1);
			
			if(b1=='>'){format=FA;}
			else if(b1=='@'){
				if(b3=='+'){format=FQ;}
				else if(b2<0 || b2=='@'){format=SAM;}
				else{format=UNKNOWN;} //probably a truncated fastq file?
			}else if(b1=='#'){
				if(s1.startsWith("#SZ:") || s1.startsWith("#SIZE:")){
					format=SKETCH;
					int x1=s1.indexOf(':');
					int x2=s1.indexOf('\t');
					if(x2>x1){
						try {
							len=Integer.parseInt(s1.substring(x1+1, x2));
						} catch (NumberFormatException e) {}
					}
				}else if(s1.startsWith("#FASTR") || s1.startsWith("#FR")){
					format=FASTR;
					if(s1.endsWith("\tINT")){i=INTERLEAVED;}
					else{i=SINGLE;}
				}else if(s1.startsWith("##fileformat=VCF")){
					format=VCF;
				}else if(s1.startsWith("#fileformat\tVar_")){
					format=VAR;
				}else if(s1.startsWith("##gff-version")){
					format=GFF;
				}else if(s1.startsWith("LOCUS ")){
					format=GBFF;
				}else if(s1.startsWith("#id\tnumericID\t")){
					format=BREAD;
				}else{format=TEXT;}
			}
//			else{format=BREAD;} //or possibly scarf

			if(format!=FQ){len=-1;}
		}
		
		if(format==FQ || forceFastq){
			boolean oldDQ=stream.FASTQ.DETECT_QUALITY;
			byte oldQin=stream.FASTQ.ASCII_OFFSET;
			byte oldQout=stream.FASTQ.ASCII_OFFSET_OUT;
			stream.FASTQ.DETECT_QUALITY=true;
			q=stream.FASTQ.testQuality(oct);
			i=(byte)(stream.FASTQ.testInterleaved(oct, fname, false) ? INTERLEAVED : SINGLE);
			//		stream.FASTQ.DETECT_QUALITY=old;
			{
				String h=oct.size()>0 ? oct.get(0) : null;
				String a=oct.size()>1 ? oct.get(1) : null;
				String b=oct.size()>5 ? oct.get(5) : null;
				if(a!=null){len=Tools.max(a.length(), len);}
				if(b!=null){len=Tools.max(b.length(), len);}
				if(len<2){len=-1;}
				if(h!=null){
					String barcode=Read.headerToBarcode(h, false);
					if(barcode!=null && barcode.length()>0) {
						numBarcodes=1;
						if(barcode.length()>2){
							barcodeDelimiter=barcodeDelimiter(barcode);
							numBarcodes=(countLetters(barcode))>0 ? (barcodeDelimiter>0 ? 2 : 1) : 0;
						}
						if(numBarcodes==1) {
							bcLen1=barcode.length();
							bcLen2=0;
						}else {
							bcLen1=barcode.indexOf(barcodeDelimiter);
							bcLen2=barcode.length()-bcLen1-1;
						}
					}
				}
			}
			stream.FASTQ.DETECT_QUALITY=oldDQ;
			stream.FASTQ.ASCII_OFFSET=oldQin;
			stream.FASTQ.ASCII_OFFSET_OUT=oldQout;
		}
		int[] r=new int[] {q, i, len, format, numBarcodes, barcodeDelimiter, bcLen1, bcLen2};
		if(verbose){System.err.println(Arrays.toString(r));}
		return r;
	}
	
	public static byte barcodeDelimiter(String barcode){
		if(barcode==null || barcode.length()<3) {return 0;}
		int letters=0, nonletters=0;
		byte delimiter=0;
		for(int i=0; i<barcode.length(); i++){
			char c=barcode.charAt(i);
			if(Tools.isLetter(c)){
				letters++;
			}else{
				if(nonletters==0) {delimiter=(byte)c;}
				nonletters++;
			}
		}
		if(nonletters==1 && letters>1){
			return delimiter;
		}
		return 0;//No delimiter or multiple delimiters
	}
	
	private static int countLetters(String s){
		int letters=0;
		for(int i=0; s!=null && i<s.length(); i++){
			char c=s.charAt(i);
			if(Tools.isLetter(c)){
				letters++;
			}
		}
		return letters;
	}
	
	public static boolean isFastaExt(String ext){
		if(ext==null){return false;}
		return (ext.equals("fa") || ext.equals("fasta") || ext.equals("fas") || ext.equals("fna") || ext.equals("ffn")
			|| ext.equals("frn") || ext.equals("seq") || ext.equals("fsa") || ext.equals("faa") || ext.equals("prot"));
	}
	
	public static boolean isFastaFile(String fname){
		if(fname==null){return false;}
		String ext=ReadWrite.rawExtension(fname);
		return isFastaExt(ext);
	}
	
	public static boolean isPgmFile(String fname){
		if(fname==null){return false;}
		String ext=ReadWrite.rawExtension(fname);
		return isPgmExt(ext);
	}
	
	public static boolean isAminoExt(String ext){
		if(ext==null){return false;}
		return ext.equals("faa") || ext.equals("prot"); //TODO: Investigate whether other extensions imply AA.
	}
	
	public static boolean isStdin(String s) {
		if(s==null){return false;}
//		if(new File(s).exists()){return false;}
		int idx=s.indexOf('.');
		if(idx>=0){s=s.substring(0, idx);}
		return s.equalsIgnoreCase("stdin") || s.equalsIgnoreCase("standardin");
	}
	public static boolean isStdout(String s) {
		if(s==null){return false;}
//		if(new File(s).exists()){return false;}
		int idx=s.indexOf('.');
		if(idx>=0){s=s.substring(0, idx);}
		return s.equalsIgnoreCase("stdout") || s.equalsIgnoreCase("standardout");
	}
	public static boolean isStderr(String s) {
		if(s==null){return false;}
//		if(new File(s).exists()){return false;}
		int idx=s.indexOf('.');
		if(idx>=0){s=s.substring(0, idx);}
		return s.equalsIgnoreCase("stderr") || s.equalsIgnoreCase("standarderror");
	}
	public static boolean isStdio(String s) {
		if(s==null){return false;}
//		if(new File(s).exists()){return false;}
		int idx=s.indexOf('.');
		if(idx>=0){s=s.substring(0, idx);}
		return s.equalsIgnoreCase("stdin") || s.equalsIgnoreCase("standardin") ||
				s.equalsIgnoreCase("stdout") || s.equalsIgnoreCase("standardout") ||
				s.equalsIgnoreCase("stderr") || s.equalsIgnoreCase("standarderror");
	}
	
//	public static boolean isStdio(String s){
//		if(s==null){return false;}
//		if(new File(s).exists()){return false;}
//		if(s.contains(".")){s=s.substring(0, s.indexOf('.'));
//		}
//		return (s.equalsIgnoreCase("stdin") || s.equalsIgnoreCase("stdout") || s.equalsIgnoreCase("stderr"));
//	}
	
	public static boolean isFastqExt(String ext){
		if(ext==null){return false;}
		return (ext.equalsIgnoreCase("fq") || ext.equalsIgnoreCase("fastq"));
	}
	
	public static boolean isPgmExt(String ext){
		if(ext==null){return false;}
		return (ext.equalsIgnoreCase("pgm") || ext.equalsIgnoreCase("pkm"));
	}
	
	public static boolean isFastqFile(String fname){
		if(fname==null){return false;}
		String ext=ReadWrite.rawExtension(fname);
		return isFastqExt(ext);
	}
	
	public static boolean isSamOrBamExt(String ext){
		if(ext==null){return false;}
		return (ext.equalsIgnoreCase("sam") || ext.equalsIgnoreCase("bam"));
	}
	
	public static boolean isGffExt(String ext){
		if(ext==null){return false;}
		return (ext.equalsIgnoreCase("gff"));
	}
	
	public static boolean isSamOrBamFile(String fname){
		if(fname==null){return false;}
		String ext=ReadWrite.rawExtension(fname);
		return isSamOrBamExt(ext);
	}
	
	public static boolean isBamExt(String ext){
		if(ext==null){return false;}
		return ext.equalsIgnoreCase("bam");
	}
	
	public static boolean isVcfExt(String ext){
		if(ext==null){return false;}
		return ext.equalsIgnoreCase("vcf");
	}
	
	public static boolean isBamFile(String fname){
		if(fname==null){return false;}
		String ext=ReadWrite.rawExtension(fname);
		return isBamExt(ext);
	}
	
	public static boolean isVcfFile(String fname){
		if(fname==null){return false;}
		String ext=ReadWrite.rawExtension(fname);
		return isVcfExt(ext);
	}

	public void deleteIfPresent() {
		File f=new File(name);
		if(f.exists()){f.delete();}
	}

	public int barcodeDelimiter() {return barcodeDelimiter;}
	public int barcodesPerRead() {return numBarcodes;}
	public int barcodeLength(int barcodeNum) {
		return barcodeNum==1 ? barcodeLength1 : barcodeNum==2 ? barcodeLength2 : 0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean extensionEquals(String ext) {
		String s=ReadWrite.rawExtension(name);
		if(s==ext) {return true;}
		if(s==null || ext==null) {return false;}
		return ext.equals(s);
	}
	
	public String rawExtension() {
		return ReadWrite.rawExtension(name);
	}
	public int rawExtensionCode() {
		String ext=ReadWrite.rawExtension(name);
		String comp=ReadWrite.compressionType(name);
		return rawExtensionCode(ext, comp);
	}
	private int rawExtensionCode(String ext, String comp) {
		if(ext==null){return UNKNOWN;}
		else if(ext.equals("fq") || ext.equals("fastq") || (comp!=null && comp.equals("fqz"))){return FASTQ;}
		else if(isFastaExt(ext)){return FASTA;}
		else if(ext.equals("bread")){return BREAD;}
		else if(ext.equals("sam")){return SAM;}
		else if(ext.equals("csfasta")){return CSFASTA;}
		else if(ext.equals("qual")){return QUAL;}
		else if(ext.equals("bam")){return BAM;}
		else if(ext.equals("sites") || ext.equals("sitesonly")){return SITES;}
		else if(ext.equals("info") || ext.equals("attachment")){return ATTACHMENT;}
		else if(ext.equals("scarf")){return SCARF;}
		else if(ext.equals("phylip")){return PHYLIP;}
		else if(ext.equals("header") || ext.equals("headers")){return HEADER;}
		else if(ext.equals("int1d")){return INT1D;}
		else if(ext.equals("long1d")){return LONG1D;}
		else if(ext.equals("bitset")){return BITSET;}
		else if(ext.equals("sketch")){return SKETCH;}
		else if(ext.equals("oneline") || ext.equals("flat")){return ONELINE;}
		else if(ext.equals("fastr") || ext.equals("fr")){return FASTR;}
		else if(ext.equals("vcf")){return VCF;}
		else if(ext.equals("var")){return VAR;}
		else if(ext.equals("gff") || ext.equals("gff3")){return GFF;}
		else if(ext.equals("bed")){return BED;}
		else if(ext.equals("pgm") || ext.equals("pkm")){return PGM;}
		else if(ext.equals("embl")){return EMBL;}
		else if(ext.equals("gbk")){return GBK;}
		else if(ext.equals("gbff")){return GBFF;}
		else if(ext.equals("bbnet")){return BBNET;}
		else if(ext.equals("bbvec") || ext.equals("vec")){return BBVEC;}
		else if(ext.equals("clade") || ext.equals("spectra")){return CLADE;}
		else if(ext.equals("txt") || ext.equals("text") || ext.equals("tsv") || ext.equals("csv")){return TXT;}
		return UNKNOWN;
	}

	public final String name(){return name;}
	public final String simpleName(){return simpleName;}
	public final int format(){return format;}
	public final int compression(){return compression;}
	public final int type(){return type;}
	public final int mode(){return mode;}
	public final boolean amino(){return amino;}
	public final boolean hasName(){return name!=null;}
	public final int asciiOffset(){return asciiOffset;}
	public final int length(){return length;}
	
	public final boolean canWrite(){
		assert(write());
		if(stdio() || devnull()){return true;}
		assert(hasName());
		File f=new File(name);
		if(!f.exists()){return true;}
		if(!f.canWrite()){return false;}
		return overwrite() || append();
	}
	
	public final boolean canRead(){
		assert(read());
		if(stdio()){return true;}
		assert(hasName());
		File f=new File(name);
		return f.canRead();
	}
	
	public final boolean unknownField(){return unknownFormat() || unknownCompression() || unknownType() || unknownMode();}

	public final boolean unknownFormat(){return format<=UNKNOWN;}
	public final boolean fasta(){return format==FASTA;}
	public final boolean fastq(){return format==FASTQ;}
	public final boolean fastr(){return format==FASTR;}
	public final boolean bread(){return format==BREAD;}
	public final boolean sam(){return format==SAM;}
	public final boolean samOrBam(){return format==SAM || format==BAM;}
	public final boolean csfasta(){return format==CSFASTA;}
	public final boolean qual(){return format==QUAL;}
	public final boolean sequential(){return format==SEQUENTIAL;}
	public final boolean random(){return format==RANDOM;}
	public final boolean sites(){return format==SITES;}
	public final boolean attachment(){return format==ATTACHMENT;}
	public final boolean header(){return format==HEADER;}
	public final boolean bam(){return format==BAM;}
	public final boolean scarf(){return format==SCARF;}
	public final boolean text(){return format==TEXT;}
	public final boolean int1d(){return format==INT1D;}
	public final boolean long1d(){return format==LONG1D;}
	public final boolean bitset(){return format==BITSET;}
	public final boolean sketch(){return format==SKETCH;}
	public final boolean oneline(){return format==ONELINE;}
	public final boolean var(){return format==VAR;}
	public final boolean vcf(){return format==VCF;}
	public final boolean gff(){return format==GFF;}
	public final boolean bed(){return format==BED;}
	public final boolean pgm(){return format==PGM;}
	public final boolean embl(){return format==EMBL;}
	public final boolean gbk(){return format==GBK;}
	public final boolean gbff(){return format==GBFF;}
	public final boolean alm(){return format==ALM;}
	public final boolean bbnet(){return format==BBNET;}
	public final boolean bbvec(){return format==BBVEC;}
	public final boolean clade(){return format==CLADE;}
	
	public final boolean preferShreds(){
		return preferShreds;
	}
	
	public boolean isSequence() {return fasta() || fastq() || fastr() || bread() || samOrBam() || csfasta() || scarf() || header() || oneline() || gbk() || embl();}

	public final boolean unknownCompression(){return compression<=UNKNOWN;}
	public final boolean raw(){return compression==RAW;}
	public final boolean gzip(){return compression==GZIP;}
	public final boolean zstd(){return compression==ZSTD;}
	public final boolean zip(){return compression==ZIP;}
	public final boolean bz2(){return compression==BZ2;}
	public final boolean fqz(){return compression==FQZ;}
	public final boolean lz(){return compression==LZ;}
	public final boolean xz(){return compression==XZ;}
	public final boolean sevenz(){return compression==SEVENZ;}
	public final boolean dsrc(){return compression==DSRC;}
	public final boolean compressed(){return compression!=RAW || format==BAM;}

	public final boolean unknownType(){return type<=UNKNOWN;}
	public final boolean file(){return type==FILE;}
	public final boolean stdio(){return type==STDIO;}
	public final boolean stdin(){return type==STDIO && input;}
	public final boolean stdout(){return type==STDIO && !input;}
	public final boolean devnull(){return type==DEVNULL;}

	public final boolean unknownMode(){return mode<=UNKNOWN;}
	public final boolean read(){return mode==READ;}
	public final boolean write(){return mode==WRITE;}

	public final boolean overwrite(){return overwrite;}
	public final boolean append(){return append;}
	public final boolean allowSubprocess(){return allowSubprocess;}
	public final boolean ordered(){return ordered;}

	public boolean interleaved(){return interleaving==INTERLEAVED;}
	
	public final boolean exists(){
		if(!file()){return read();}
		File f=new File(name);
		if(!f.exists() && !gzip()){return false;}
		long size=f.length();
		return size>10;
	}
	
//	public final boolean interleaved(){return interleaved;}
//	public final long maxReads(){return maxReads;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String name;
	private final String simpleName;
	private final int format;
	private final int asciiOffset;
	private final int compression;
	private final int type;
	private final int mode;
	private final int interleaving;
	private final int length;
	private final boolean input;
	private final boolean amino;
	private final int numBarcodes;
	private final int barcodeDelimiter;
	private final int barcodeLength1;
	private final int barcodeLength2;

	private final boolean overwrite;
	private final boolean append;
	private final boolean allowSubprocess;
	private final boolean ordered;
	
//	private final int magicNumber;
	
	public boolean preferShreds=false;
//	private final long maxReads;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean verbose=false;
	public static boolean PRINT_WARNING=true;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	public static final int UNKNOWN=0;
	
	/* Format */
	
	public static final int FA=1, FASTA=1;
	public static final int FQ=2, FASTQ=2;
	public static final int BREAD=3;
	public static final int SAM=4;
	public static final int CSFASTA=5;
	public static final int QUAL=6;
	public static final int SEQUENTIAL=7;
	public static final int RANDOM=8;
	public static final int SITES=9;
	public static final int ATTACHMENT=10;
	public static final int BAM=11;
	public static final int SCARF=12;
	public static final int TEXT=13, TXT=13;
	public static final int PHYLIP=14;
	public static final int HEADER=15;
	public static final int INT1D=16;
	public static final int LONG1D=17;
	public static final int BITSET=18;
	public static final int SKETCH=19;
	public static final int ONELINE=20;
	public static final int FR=21, FASTR=21;
	public static final int VCF=22;
	public static final int VAR=23;
	public static final int GFF=24;
	public static final int BED=25;
	public static final int PGM=26, PKM=26;
	public static final int EMBL=27;
	public static final int GBK=28;
	public static final int GBFF=29;//TODO: this may be the same as GBK...
	//Alignment Model, from Consensus package
	public static final int ALM=30;
	public static final int BBNET=31;
	public static final int BBVEC=32;
	public static final int CLADE=33;
	
	public static final String[] FORMAT_ARRAY=new String[] {
		"unknown", "fasta", "fastq", "bread", "sam", "csfasta",
		"qual", "sequential", "random", "sites", "attachment",
		"bam", "scarf", "text", "phylip", "header", "int1d",
		"long1d", "bitset", "sketch", "oneline", "fastr",
		"vcf", "var", "gff", "bed", "pgm", "embl", "gbk", "gbff",
		"alm", "bbnet", "bbvec", "vec", "clade", "spectra"
	};
	
	public static final String[] EXTENSION_LIST=new String[] {
		"fq", "fastq", "fa", "fasta", "fas", "fna",
		"ffn", "frn", "seq", "fsa", "faa",
		"bread", "sam", "csfasta", "qual", "bam",
		"scarf", "phylip", "txt",
		"gz", "gzip", "bz2", "zip", "xz", "dsrc", "header", "headers",
		"int1d", "long1d", "bitset", "sketch", "oneline", "flat", "fqz",
		"gff", "gff3", "var", "vcf", "bed", "pgm", "embl", "gbk", "gbff", "alm", 
		"bbnet", "bbvec", "vec", "clade", "spectra", "7z", "zst"
	};
	
	/* Compression */
	
	public static final int RAW=1;
	public static final int GZ=2, GZIP=2;
	public static final int ZIP=3;
	public static final int BZ2=4;
	public static final int XZ=5;
	public static final int c4=6;
	public static final int SEVENZ=7;
	public static final int DSRC=8;
	public static final int FQZ=9;
	public static final int LZ=10;
	public static final int AC=11;
	public static final int ZSTD=12;
	
	public static final String[] COMPRESSION_ARRAY=new String[] {
		"unknown", "raw", "gz", "zip", "bz2", "xz",
		"c4", "7z", "dsrc", "fqz", "lz", "ac", "zst"
	};
	
	/* Type */
	
	public static final int FILE=1;
	public static final int STDIO=2, STDIN=2, STDOUT=2;
	public static final int DEVNULL=3;
//	public static final int NULL=4;
	
	private static final String[] TYPE_ARRAY=new String[] {
		"unknown", "file", "stdio", "devnull"
	};
	
	/* Mode */
	
	public static final int READ=1, WRITE=2;
	
	private static final String[] MODE_ARRAY=new String[] {
		"unknown", "read", "write"
	};
	
	/* Interleaving */
	
	public static final int SINGLE=1, INTERLEAVED=2;
	
	private static final String[] INTERLEAVING_ARRAY=new String[] {
		"unknown", "single-ended", "interleaved"
	};
	
}
