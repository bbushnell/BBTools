package driver;

import java.io.OutputStream;
import java.io.PrintWriter;

import fileIO.ReadWrite;
import fileIO.TextFile;

/**
 * Converts SAM alignment files to simplified ALN format.
 * Extracts chromosome, position, and orientation from SAM records.
 * Filters out unmapped reads and compresses output with gzip.
 * @author Brian Bushnell
 */
public class ConvertSamToAln {
	
	/**
	 * Program entry point that converts multiple SAM files.
	 * Processes each filename argument and reports completion status.
	 * @param args SAM filenames to convert
	 */
	public static void main(String[] args){
		for(String s : args){
			convert(s);
			System.out.println("Converted "+s);
		}
	}
	
	/**
	 * Converts a single SAM file to ALN format.
	 * Reads SAM records, extracts alignment information for mapped reads,
	 * and writes chromosome, position, and strand to compressed ALN file.
	 * Output filename replaces .sam extension with .aln.gz.
	 *
	 * @param fname Input SAM filename (may be compressed)
	 */
	public static final void convert(String fname){
		TextFile tf=new TextFile(fname, false);
		
		
		
		String outname=fname;
		if(outname.toLowerCase().endsWith(".zip")){outname=outname.substring(0, outname.length()-4);}
		if(outname.toLowerCase().endsWith(".gz")){outname=outname.substring(0, outname.length()-3);}
		if(outname.toLowerCase().endsWith(".bz2")){outname=outname.substring(0, outname.length()-4);}
		if(outname.toLowerCase().endsWith(".sam")){outname=outname.substring(0, outname.length()-4);}
		outname=outname+".aln.gz";
		
		String s=null;
		
		OutputStream os=ReadWrite.getOutputStream(outname, false, true, true);
		PrintWriter out=new PrintWriter(os);
		
		for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(!s.startsWith("@")){
				String[] line=s.split("\t");
				assert(line.length>1) : s;
				
				boolean success=true;
				boolean nomap=false;
				boolean reverse=false;
				
				int flag=-1;
				String chrom=null;
				int loc=-1;
				
				try {
					flag=Integer.parseInt(line[1]);
					chrom=line[2];
					loc=Integer.parseInt(line[3]);
					nomap=((flag&0x4)!=0);
					reverse=((flag&0x10)!=0);
				} catch (NumberFormatException e) {
					success=false;
				}
				
				if(success && !nomap){
					String aln=chrom+"\t"+loc+"\t"+(reverse ? "R" : "F")+"\n";
					out.print(aln);
				}
				
				
			}
		}

		tf.close();
		out.flush();
		out.close();
		
	}
	
}
