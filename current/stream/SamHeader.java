package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import dna.Data;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Builds SAM header text used by SAM and BAM writers from shared reference/configuration data.
 * Produces @HD, optional @SQ and @RG, and @PG records. Does not validate or escape
 * configured metadata fields. Data scaffold arrays and SamLine/Shared settings must
 * be initialized and remain stable while building; this class does not snapshot them.
 * header0B/header2B append to caller-owned storage without clearing it. Their
 * StringBuilder counterparts allocate results; print methods leave output lifecycle to the caller.
 * Newline, sorting and residual-buffer rules differ between methods as documented below.
 *
 * @author Brian Bushnell
 * @date Jul 7, 2014
 */
public class SamHeader{

	/*--------------------------------------------------------------*/
	/*----------------        Header Builders       ----------------*/
	/*--------------------------------------------------------------*/

	/** Appends @HD with version 1.3 below SamLine.VERSION 1.4, otherwise 1.4, and SO:unsorted.
	 * @param bb Nonnull destination, not cleared; no terminal newline is appended
	 * @return The supplied builder */
	public static ByteBuilder header0B(ByteBuilder bb){
		//		if(MAKE_TOPHAT_TAGS){
		//			return new ByteBuilder("@HD\tVN:"+(VERSION<1.4f ? "1.0" : "1.4")+"\tSO:unsorted");
		//		}
		bb.append("@HD\tVN:");
		bb.append((SamLine.VERSION<1.4f ? "1.3" : "1.4"));
		bb.append("\tSO:unsorted");
		return bb;
	}

	/** Allocates @HD text without a final newline, using the same version rule as header0B.
	 * @return New builder containing the version and unsorted-order record */
	//Historical review found no external callers; retained StringBuilder counterpart of header0B.
	public static StringBuilder header0(){
		//		if(MAKE_TOPHAT_TAGS){
		//			return new StringBuilder("@HD\tVN:"+(SamLine.VERSION<1.4f ? "1.0" : "1.4")+"\tSO:unsorted");
		//		}
		StringBuilder sb=new StringBuilder("@HD\tVN:"+(SamLine.VERSION<1.4f ? "1.3" : "1.4")+"\tSO:unsorted");
		return sb;
	}

	/** Builds scaffold dictionary strings over an inclusive, upper-clamped chromosome range.
	 * Name processing uses scaffoldPrefixes/TRIM_RNAME; each length is capped at Integer.MAX_VALUE.
	 * @param minChrom First valid Data chromosome index; no lower-bound normalization
	 * @param maxChrom Last chromosome, capped at Data.numChroms
	 * @param sort Sort complete @SQ strings across this call's selected range
	 * @param newline Append LF to each returned string
	 * @return New list of newly allocated dictionary strings */
	static ArrayList<String> scaffolds(int minChrom, int maxChrom, boolean sort, boolean newline){
		final ArrayList<String> list=new ArrayList<String>(4000);
		final StringBuilder sb=new StringBuilder(1000);
		for(int i=minChrom; i<=maxChrom && i<=Data.numChroms; i++){
			final byte[][] inames=Data.scaffoldNames[i];
			for(int j=0; j<Data.chromScaffolds[i]; j++){
				final byte[] scn=inames[j];
				sb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scn==null){
					assert(false) : "scaffoldName["+i+"]["+j+"] = null";
					sb.append("null");
				}else{
					appendScafName(sb, scn);
				}
				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j])));
				//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
				//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"b"+Data.GENOME_BUILD).replace('\t', ' '));

				if(newline){sb.append('\n');}
				list.add(sb.toString());
				sb.setLength(0);
			}
		}
		if(sort){Shared.sort(list);}
		return list;
	}

	/** Allocates newline-terminated @SQ records, optionally sorted across the selected range.
	 * @param minChrom First valid Data chromosome index, inclusive
	 * @param maxChrom Last chromosome, inclusive and capped at Data.numChroms
	 * @return New builder; empty if the selected range has no scaffolds */
	public static StringBuilder header1(int minChrom, int maxChrom){
		StringBuilder sb=new StringBuilder(20000);
		if(SamLine.SORT_SCAFFOLDS){
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true, true);
			for(int i=0; i<scaffolds.size(); i++){
				sb.append(scaffolds.get(i));
				scaffolds.set(i, null);
			}
			return sb;
		}

		for(int i=minChrom; i<=maxChrom && i<=Data.numChroms; i++){
			final byte[][] inames=Data.scaffoldNames[i];
			for(int j=0; j<Data.chromScaffolds[i]; j++){
				byte[] scn=inames[j];
				sb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scn==null){
					assert(false) : "scaffoldName["+i+"]["+j+"] = null";
					sb.append("null");
				}else{
					appendScafName(sb, scn);
				}

				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j])));
				//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
				//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"build "+Data.GENOME_BUILD).replace('\t', ' '));

				sb.append('\n');
			}
		}

		return sb;
	}

	/** Prints newline-terminated @SQ records without flushing, closing or checking the writer.
	 * @param minChrom First valid Data chromosome index, inclusive
	 * @param maxChrom Last chromosome, inclusive and capped at Data.numChroms
	 * @param pw Borrowed nonnull destination */
	public static void printHeader1(int minChrom, int maxChrom, PrintWriter pw){
		if(SamLine.SORT_SCAFFOLDS){
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true, true);
			for(int i=0; i<scaffolds.size(); i++){
				pw.print(scaffolds.set(i, null));
			}
			return;
		}

		for(int i=minChrom; i<=maxChrom && i<=Data.numChroms; i++){
			final byte[][] inames=Data.scaffoldNames[i];
			StringBuilder sb=new StringBuilder(256);
			for(int j=0; j<Data.chromScaffolds[i]; j++){
				final byte[] scn=inames[j];
				//				StringBuilder sb=new StringBuilder(7+(scn==null ? 4 : scn.length)+4+10+4+/*(Data.name==null ? 0 : Data.name.length()+1)+11*/+4);//last one could be 1
				sb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scn==null){
					assert(false) : "scaffoldName["+i+"]["+j+"] = null";
					sb.append("null");
				}else{
					appendScafName(sb, scn);
				}
				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j])));
				//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
				//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"b"+Data.GENOME_BUILD).replace('\t', ' '));

				sb.append('\n');

				pw.print(sb);
				sb.setLength(0);
			}
		}
	}

	/** Appends newline-terminated @SQ records and writes whenever bb reaches 32768 bytes.
	 * Existing bytes in bb are included in threshold writes; successful writes clear it.
	 * The final residue stays in bb for the caller to append/flush. Does not flush or
	 * close os. SORT_SCAFFOLDS sorts across this call's range; sorted-path I/O errors
	 * use KillSwitch.exceptionKill, while the unsorted path wraps them in RuntimeException.
	 * @param minChrom First valid Data chromosome index, inclusive
	 * @param maxChrom Last chromosome, inclusive and capped at Data.numChroms
	 * @param bb Borrowed nonnull buffer; final unwritten bytes remain here
	 * @param os Borrowed nonnull output stream for threshold writes */
	public static void printHeader1B(int minChrom, int maxChrom, ByteBuilder bb, OutputStream os){
		if(verbose){System.err.println("printHeader1B("+minChrom+", "+maxChrom+")");}

		if(SamLine.SORT_SCAFFOLDS){
			if(verbose){System.err.println("Sorting scaffolds");}
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true, true);
			for(int i=0; i<scaffolds.size(); i++){
				String s=scaffolds.set(i, null);
				bb.append(s);
				if(bb.length>=32768){
					try{
						os.write(bb.array, 0, bb.length);
					}catch(IOException e){
						KillSwitch.exceptionKill(e);
					}
					bb.setLength(0);
				}
			}
			return;
		}

		if(verbose){System.err.println("Iterating over chroms");}
		for(int chrom=minChrom; chrom<=maxChrom && chrom<=Data.numChroms; chrom++){
			//			if(verbose){System.err.println("chrom "+chrom);}
			final byte[][] inames=Data.scaffoldNames[chrom];
			//			if(verbose){System.err.println("inames"+(inames==null ? " = null" : ".length = "+inames.length));}
			final int numScafs=Data.chromScaffolds[chrom];
			//			if(verbose){System.err.println("scaffolds: "+numScafs);}
			assert(inames.length==numScafs) : "Mismatch between number of scaffolds and names for chrom "+chrom+": "+inames.length+" != "+numScafs;
			for(int scaf=0; scaf<numScafs; scaf++){
				//				if(verbose){System.err.println("chromScaffolds["+scaf+"] = "+(inames==null ? "=null" : ".length="+inames.length));}
				final byte[] scafName=inames[scaf];
				//				if(verbose){System.err.println("scafName = "+(scafName==null ? "null" : new String(scafName)));}
				bb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scafName==null){
					assert(false) : "scaffoldName["+chrom+"]["+scaf+"] = null";
					//[stream/SamHeader#002 corrected] With assertions off, ByteBuilder.append(byte[])
					//writes literal "null", matching the StringBuilder paths; the old NPE claim was stale.
					bb.append(scafName);
				}else{
					appendScafName(bb, scafName);
				}
				bb.append("\tLN:");
				bb.append(Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[chrom][scaf])));
				//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
				//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"b"+Data.GENOME_BUILD).replace('\t', ' '));

				bb.nl();

				if(bb.length>=32768){
					try{
						os.write(bb.array, 0, bb.length);
					}catch(IOException e){
						throw new RuntimeException(e);
					}
					bb.setLength(0);
				}
			}
		}
	}

	/** Submits newline-terminated @SQ text to an already started TextStreamWriter.
	 * The sorted path submits immutable strings; the unsorted path has the mutable
	 * buffer ownership issue annotated below. Does not finish or close the writer.
	 * @param minChrom First valid Data chromosome index, inclusive
	 * @param maxChrom Last chromosome, inclusive and capped at Data.numChroms
	 * @param tsw Borrowed started text writer */
	public static void printHeader1(int minChrom, int maxChrom, TextStreamWriter tsw){
		if(SamLine.SORT_SCAFFOLDS){
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true, true);
			for(int i=0; i<scaffolds.size(); i++){
				tsw.print(scaffolds.set(i, null));
			}
			return;
		}

		for(int i=minChrom; i<=maxChrom && i<=Data.numChroms; i++){
			final byte[][] inames=Data.scaffoldNames[i];
			final StringBuilder sb=new StringBuilder(256);
			for(int j=0; j<Data.chromScaffolds[i]; j++){
				final byte[] scn=inames[j];
				//				StringBuilder sb=new StringBuilder(7+(scn==null ? 4 : scn.length)+4+10+4+/*(Data.name==null ? 0 : Data.name.length()+1)+11*/+4);//last one could be 1
				sb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scn==null){
					assert(false) : "scaffoldName["+i+"]["+j+"] = null";
					sb.append("null");
				}else{
					appendScafName(sb, scn);
				}
				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j])));
				//				sb.append("\tLN:"+Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[i][j]+1000L)));
				//				sb.append("\tAS:"+((Data.name==null ? "" : Data.name+" ")+"b"+Data.GENOME_BUILD).replace('\t', ' '));

				sb.append('\n');

				//TODO: Probable bug - TextStreamWriter.print(CharSequence) retains this builder,
				//then setLength(0) below clears queued text before asynchronous consumption.
				//No callers of this legacy overload found in the current source tree; fix separately.
				tsw.print(sb);
				sb.setLength(0);
			}
		}
	}

	/** Appends a processed scaffold name by casting retained bytes to chars.
	 * When scaffoldPrefixes is enabled, skips through the first dollar sign; absence
	 * asserts when enabled and otherwise yields an empty name. TRIM_RNAME stops at
	 * the first Character.isWhitespace byte. Does not clear or terminate the builder.
	 * @param sb Nonnull destination
	 * @param scn Nonnull scaffold-name bytes */
	static void appendScafName(StringBuilder sb, byte[] scn){
		int start=0;
		if(Data.scaffoldPrefixes){
			while(start<scn.length && scn[start]!='$'){start++;}
			start++;
			assert(start<=scn.length) : "Scaffold name missing '$' prefix delimiter: "+new String(scn);
			if(start>scn.length){start=scn.length;}
		}
		int end=scn.length;
		if(Shared.TRIM_RNAME){
			for(int i=start; i<end; i++){
				if(Character.isWhitespace(scn[i])){end=i; break;}
			}
		}
		final int len=end-start;
		if(len>0){
			final char[] buffer=Shared.getTLCB(len);
			for(int i=start; i<end; i++){buffer[i-start]=(char)scn[i];}
			sb.append(buffer, 0, len);
		}
	}

	/** Appends the retained scaffold-name bytes with the StringBuilder overload's bounds rules.
	 * Prefix stripping and whitespace trimming select a byte slice without decoding it.
	 * @param sb Nonnull destination, not cleared or newline-terminated
	 * @param scn Nonnull scaffold-name bytes */
	static void appendScafName(ByteBuilder sb, byte[] scn){
		int start=0;
		if(Data.scaffoldPrefixes){
			while(start<scn.length && scn[start]!='$'){start++;}
			start++;
			assert(start<=scn.length) : "Scaffold name missing '$' prefix delimiter: "+new String(scn);
			if(start>scn.length){start=scn.length;}
		}
		int end=scn.length;
		if(Shared.TRIM_RNAME){
			for(int i=start; i<end; i++){
				if(Character.isWhitespace(scn[i])){end=i; break;}
			}
		}
		final int len=end-start;
		if(len>0){sb.append(scn, start, len);}
	}

	/** Allocates optional @RG plus @PG text from current global settings.
	 * The @RG record, when enabled by READGROUP_ID, ends with LF; @PG has no final newline.
	 * Command text includes JVM/command arguments verbatim and always prefixes a
	 * nonnull BBMAP_CLASS with " align2.", unlike header2B's leading-space handling.
	 * @return New builder containing read-group/program metadata */
	public static StringBuilder header2(){
		StringBuilder sb=new StringBuilder(1000);
		//		sb.append("@RG\tID:unknownRG\tSM:unknownSM\tPL:ILLUMINA\n"); //Can cause problems.  If RG is in the header, reads may need extra fields.

		//		if(MAKE_TOPHAT_TAGS){
		////			sb.append("@PG\tID:TopHat\tVN:2.0.6\tCL:/usr/common/jgi/aligners/tophat/2.0.6/bin/tophat -p 16 -r 0 --max-multihits 1 Creinhardtii_236 reads_1.fa reads_2.fa");
		//			sb.append("@PG\tID:TopHat\tVN:2.0.6");
		//		}else{
		//			sb.append("@PG\tID:BBMap\tPN:BBMap\tVN:"+Shared.BBMAP_VERSION_STRING);
		//		}

		if(SamLine.READGROUP_ID!=null){
			sb.append("@RG\tID:").append(SamLine.READGROUP_ID);
			if(SamLine.READGROUP_CN!=null){sb.append("\tCN:").append(SamLine.READGROUP_CN);}
			if(SamLine.READGROUP_DS!=null){sb.append("\tDS:").append(SamLine.READGROUP_DS);}
			if(SamLine.READGROUP_DT!=null){sb.append("\tDT:").append(SamLine.READGROUP_DT);}
			if(SamLine.READGROUP_FO!=null){sb.append("\tFO:").append(SamLine.READGROUP_FO);}
			if(SamLine.READGROUP_KS!=null){sb.append("\tKS:").append(SamLine.READGROUP_KS);}
			if(SamLine.READGROUP_LB!=null){sb.append("\tLB:").append(SamLine.READGROUP_LB);}
			if(SamLine.READGROUP_PG!=null){sb.append("\tPG:").append(SamLine.READGROUP_PG);}
			if(SamLine.READGROUP_PI!=null){sb.append("\tPI:").append(SamLine.READGROUP_PI);}
			if(SamLine.READGROUP_PL!=null){sb.append("\tPL:").append(SamLine.READGROUP_PL);}
			if(SamLine.READGROUP_PU!=null){sb.append("\tPU:").append(SamLine.READGROUP_PU);}
			if(SamLine.READGROUP_SM!=null){sb.append("\tSM:").append(SamLine.READGROUP_SM);}
			sb.append('\n');
		}

		sb.append("@PG\tID:BBMap\tPN:"+PN+"\tVN:");
		sb.append(Shared.BBTOOLS_VERSION_STRING);
//		assert(false) : sb+"\n"+PN;
		if(Shared.BBMAP_CLASS!=null){
			sb.append("\tCL:java");
			{
				List<String> list=null;
				list=Shared.JVM_ARGS();
				if(list!=null){
					for(String s : list){
						sb.append(' ');
						sb.append(s);
					}
				}
			}
			sb.append(" align2."+Shared.BBMAP_CLASS);
			if(Shared.COMMAND_LINE!=null){
				for(String s : Shared.COMMAND_LINE){
					sb.append(' ');
					sb.append(s);
				}
			}
		}

		return sb;
	}

	/** Appends optional @RG plus @PG text without clearing the destination.
	 * The @RG record ends with LF; @PG does not. Uses ID:BBMap and configurable PN. Command text
	 * uses current JVM/command arguments verbatim; a BBMAP_CLASS beginning with a space
	 * suppresses the usual " align2." prefix. No metadata escaping is performed.
	 * @param sb Nonnull borrowed destination
	 * @return The supplied builder */
	public static ByteBuilder header2B(ByteBuilder sb){

		if(SamLine.READGROUP_ID!=null){
			sb.append("@RG\tID:").append(SamLine.READGROUP_ID);
			if(SamLine.READGROUP_CN!=null){sb.append("\tCN:").append(SamLine.READGROUP_CN);}
			if(SamLine.READGROUP_DS!=null){sb.append("\tDS:").append(SamLine.READGROUP_DS);}
			if(SamLine.READGROUP_DT!=null){sb.append("\tDT:").append(SamLine.READGROUP_DT);}
			if(SamLine.READGROUP_FO!=null){sb.append("\tFO:").append(SamLine.READGROUP_FO);}
			if(SamLine.READGROUP_KS!=null){sb.append("\tKS:").append(SamLine.READGROUP_KS);}
			if(SamLine.READGROUP_LB!=null){sb.append("\tLB:").append(SamLine.READGROUP_LB);}
			if(SamLine.READGROUP_PG!=null){sb.append("\tPG:").append(SamLine.READGROUP_PG);}
			if(SamLine.READGROUP_PI!=null){sb.append("\tPI:").append(SamLine.READGROUP_PI);}
			if(SamLine.READGROUP_PL!=null){sb.append("\tPL:").append(SamLine.READGROUP_PL);}
			if(SamLine.READGROUP_PU!=null){sb.append("\tPU:").append(SamLine.READGROUP_PU);}
			if(SamLine.READGROUP_SM!=null){sb.append("\tSM:").append(SamLine.READGROUP_SM);}
			sb.append('\n');
		}

		//[stream/SamHeader#001 FIXED 2026-06-20 (greenlit)] @PG ID unified to BBMap: this (ByteBuilder/BAM path via
		//makeHeaderList) now emits ID:BBMap, matching the StringBuilder twin header2()(L289, SAM-text). Brian: SAM/BAM
		//usually originate from BBMap (or FilterSam etc.), so BBMap is the right @PG ID; the old ID:BBTools here diverged.
		sb.append("@PG\tID:BBMap\tPN:"+PN+"\tVN:");
		sb.append(Shared.BBTOOLS_VERSION_STRING);
//		assert(false) : sb+"\n"+PN;
		if(Shared.BBMAP_CLASS!=null){
			sb.append("\tCL:java");
			{
				List<String> list=null;
				list=Shared.JVM_ARGS();
				if(list!=null){
					for(String s : list){
						sb.append(' ');
						sb.append(s);
					}
				}
			}
			if(!Shared.BBMAP_CLASS.startsWith(" ")){sb.append(" align2.");}
			sb.append(Shared.BBMAP_CLASS);
			if(Shared.COMMAND_LINE!=null){
				for(String s : Shared.COMMAND_LINE){
					sb.append(' ');
					sb.append(s);
				}
			}
		}

		return sb;
	}
	/** Allocates @SQ byte arrays without terminal LF over an inclusive chromosome range.
	 * Sorted output converts scaffold strings with the platform-default charset;
	 * unsorted output copies ByteBuilder content directly. Returned arrays are independent.
	 * @param minChrom First valid Data chromosome index
	 * @param maxChrom Last chromosome, capped at Data.numChroms
	 * @return New list of newly allocated dictionary records */
	private static ArrayList<byte[]> header1B(int minChrom, int maxChrom){
		if(verbose){System.err.println("printHeader1B("+minChrom+", "+maxChrom+")");}

		if(SamLine.SORT_SCAFFOLDS){
			if(verbose){System.err.println("Sorting scaffolds");}
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true, false);
			ArrayList<byte[]> list=new ArrayList<byte[]>(scaffolds.size());
			for(int i=0; i<scaffolds.size(); i++){
				String s=scaffolds.set(i, null);
				list.add(s.getBytes());
			}
			return list;
		}

		if(verbose){System.err.println("Iterating over chroms");}
		ByteBuilder bb=new ByteBuilder();
		ArrayList<byte[]> list=new ArrayList<byte[]>();
		for(int chrom=minChrom; chrom<=maxChrom && chrom<=Data.numChroms; chrom++){
			final byte[][] inames=Data.scaffoldNames[chrom];
			final int numScafs=Data.chromScaffolds[chrom];
			assert(inames.length==numScafs) : "Mismatch between number of scaffolds and names for chrom "+chrom+": "+inames.length+" != "+numScafs;
			for(int scaf=0; scaf<numScafs; scaf++){
				final byte[] scafName=inames[scaf];
				bb.clear();
				bb.append("@SQ\tSN:");//+Data.scaffoldNames[i][j]);
				if(scafName==null){
					assert(false) : "scaffoldName["+chrom+"]["+scaf+"] = null";
					//[stream/SamHeader#002 corrected] Current ByteBuilder.append(null) emits literal
					//"null" without assertions, not the NPE claimed by the earlier annotation.
					bb.append(scafName);
				}else{
					appendScafName(bb, scafName);
				}
				bb.append("\tLN:");
				bb.append(Tools.min(Integer.MAX_VALUE, (Data.scaffoldLengths[chrom][scaf])));
				list.add(bb.toBytes());
			}
		}
		return list;
	}

	//comprehension: THE live BAM/SAM header builder (callers: BamWriter/SamWriter/SamWriterST/SamWriterST2). One byte[] per line:
	//@HD (header0B), then @SQ lines per chrom (header1B) unless supressed, then @RG/@PG (header2B split on '\n'). The shared bb is
	//reused safely: toBytes()+clear() after @HD, untouched during the @SQ loop (header1B uses its own bb), clean again for header2B.
	/** Builds an independent list containing @HD, optional @SQ, and split @RG/@PG text.
	 * Does not append record-ending newlines to returned arrays. SQ generation is called
	 * once per chromosome, so SORT_SCAFFOLDS sorts within each chromosome here.
	 * @param supressHeaderSequences Omit only the reference dictionary when true
	 * @param MINCHROM First chromosome, or -1 for one
	 * @param MAXCHROM Last chromosome, or -1 for Data.numChroms
	 * @return New list with independent byte arrays in generated header order */
	public static ArrayList<byte[]> makeHeaderList(boolean supressHeaderSequences, int MINCHROM, int MAXCHROM){

		ArrayList<byte[]> list=new ArrayList<byte[]>(32);

		ByteBuilder bb=new ByteBuilder(4096);
		header0B(bb);
		list.add(bb.toBytes());
		bb.clear();
		int a=(MINCHROM==-1 ? 1 : MINCHROM);
		int b=(MAXCHROM==-1 ? Data.numChroms : MAXCHROM);
		if(!supressHeaderSequences){
			for(int chrom=a; chrom<=b; chrom++){
				ArrayList<byte[]> list2=header1B(chrom, chrom);
				list.addAll(list2);
			}
		}
		header2B(bb);
		byte[][] h2b=bb.split('\n');
		for(byte[] line : h2b){list.add(line);}
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Program-name value inserted into @PG; configure before building headers. */
	public static String PN="BBMap";
	/** Enables diagnostic header-generation logging when compiled true. */
	private static final boolean verbose=false;

}
