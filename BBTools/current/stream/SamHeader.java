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
 * @author Brian Bushnell
 * @date Jul 7, 2014
 *
 */
public class SamHeader {

	public static ByteBuilder header0B(ByteBuilder bb){
		//		if(MAKE_TOPHAT_TAGS){
		//			return new ByteBuilder("@HD\tVN:"+(VERSION<1.4f ? "1.0" : "1.4")+"\tSO:unsorted");
		//		}
		bb.append("@HD\tVN:");
		bb.append((SamLine.VERSION<1.4f ? "1.3" : "1.4"));
		bb.append("\tSO:unsorted");
		return bb;
	}

	public static StringBuilder header0(){
		//		if(MAKE_TOPHAT_TAGS){
		//			return new StringBuilder("@HD\tVN:"+(SamLine.VERSION<1.4f ? "1.0" : "1.4")+"\tSO:unsorted");
		//		}
		StringBuilder sb=new StringBuilder("@HD\tVN:"+(SamLine.VERSION<1.4f ? "1.3" : "1.4")+"\tSO:unsorted");
		return sb;
	}

	static ArrayList<String> scaffolds(int minChrom, int maxChrom, boolean sort){
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

				sb.append('\n');
				list.add(sb.toString());
				sb.setLength(0);
			}
		}
		if(sort){Shared.sort(list);}
		return list;
	}

	public static StringBuilder header1(int minChrom, int maxChrom){
		StringBuilder sb=new StringBuilder(20000);
		if(SamLine.SORT_SCAFFOLDS){
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true);
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

	public static void printHeader1(int minChrom, int maxChrom, PrintWriter pw){
		if(SamLine.SORT_SCAFFOLDS){
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true);
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

	public static void printHeader1B(int minChrom, int maxChrom, ByteBuilder bb, OutputStream os){
		if(verbose){System.err.println("printHeader1B("+minChrom+", "+maxChrom+")");}

		if(SamLine.SORT_SCAFFOLDS){
			if(verbose){System.err.println("Sorting scaffolds");}
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true);
			for(int i=0; i<scaffolds.size(); i++){
				String s=scaffolds.set(i, null);
				bb.append(s);
				if(bb.length>=32768){
					try {
						os.write(bb.array, 0, bb.length);
					} catch (IOException e) {
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
					try {
						os.write(bb.array, 0, bb.length);
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
					bb.setLength(0);
				}
			}
		}
	}

	public static void printHeader1(int minChrom, int maxChrom, TextStreamWriter tsw){
		if(SamLine.SORT_SCAFFOLDS){
			ArrayList<String> scaffolds=scaffolds(minChrom, maxChrom, true);
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

				tsw.print(sb);
				sb.setLength(0);
			}
		}
	}

	static void appendScafName(StringBuilder sb, byte[] scn){
		if(Data.scaffoldPrefixes){
			int k=0;
			while(k<scn.length && scn[k]!='$'){k++;}
			k++;
			while(k<scn.length){
				sb.append((char)scn[k]);
				k++;
			}
		}else{
			final char[] buffer=Shared.getTLCB(scn.length);
			for(int i=0; i<scn.length; i++){buffer[i]=(char)scn[i];}
			sb.append(buffer, 0, scn.length);
		}
	}

	static void appendScafName(ByteBuilder sb, byte[] scn){
		if(Data.scaffoldPrefixes){
			int k=0;
			while(k<scn.length && scn[k]!='$'){k++;}
			k++;
			while(k<scn.length){
				sb.append(scn[k]);
				k++;
			}
		}else{
			sb.append(scn);
		}
	}

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
	
	public static String PN="BBMap";
	private static final boolean verbose=false;

}
