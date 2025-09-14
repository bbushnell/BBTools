package stream;

import dna.Data;
import dna.Gene;
import shared.Shared;
import structures.ByteBuilder;

/** For custom headers of synthetic reads.
 * This class is pretty ancient. */
public class CustomHeader {
	
	public CustomHeader(final String original){
		this(original, getPairnum(original));
	}
	
	public CustomHeader(final String original, int rnum_){
		assert(rnum_==0 || rnum_==1);
		rnum=rnum_;
		
		try {
			assert(original.startsWith("SYN"));
			final int middle=original.indexOf(MIDDLE);
			final int space=original.indexOf(' ');
			//if(space<0){space=original.length();}
			final String line;
			
			if(middle<0){
				line=original.substring(4, space<0 ? original.length() : space);
			}else{
				if(rnum==0){
					line=original.substring(4, middle);
				}else{
					line=original.substring(middle+1, space<0 ? original.length() : space);
				}
			}
			
			String[] split=line.split("_");
			//start_stop_insert_strand_chrom_cigar_bbstart_bbstop_bbchrom_match_rname_
			assert(split.length==9) : split.length+"\n"+line;
			
			id=Long.parseLong(split[0]);
			start=Integer.parseInt(split[1]);
			stop=Integer.parseInt(split[2]);
			insert=Integer.parseInt(split[3]);
			strand=Gene.toStrand(split[4]);
			bbstart=Integer.parseInt(split[5]);
			bbchrom=Integer.parseInt(split[6]);
			match=(split[7].equals(".") ? null : split[7].getBytes());
			rname=decodeRname(split[8]);
		} catch (Exception e) {
			FASTQ.PARSE_CUSTOM=false;
			if(FASTQ.PARSE_CUSTOM_WARNING){
				e.printStackTrace();
				System.err.println("Turned off PARSE_CUSTOM for input "+original);
			}
		}
//		System.err.println("'"+rnum+"'");
//		System.err.println("'"+bbstart+"'");
//		System.err.println(original);
//		new Exception().printStackTrace();
	}
	
	public CustomHeader(int rnum_, long id_, int start_, int stop_, int insert_, int strand_, 
			String cigar_, int bbstart_, int bbstop_, int bbchrom_, int bbscaffold_, byte[] match_, String rname_){
		assert(rnum_==0 || rnum_==1);
		rnum=rnum_;
		id=id_;
		start=start_;
		stop=stop_;
		insert=insert_;
		strand=strand_;
		bbstart=bbstart_;
		bbchrom=bbchrom_;
		match=match_;
		rname=rname_;
		
//		String s=toString();
//		Header h=new Header("SYN_"+s, rnum);
//		String s2=h.toString();
//		assert(rnum==h.rnum) : s+"\n"+s2;
//		assert(id==h.id) : s+"\n"+s2;
//		assert(start==h.start) : s+"\n"+s2;
//		assert(stop==h.stop) : s+"\n"+s2;
//		assert(insert==h.insert) : s+"\n"+s2;
//		assert(strand==h.strand) : s+"\n"+s2;
//		
//		assert(bbstart==h.bbstart) : s+"\n"+s2;
//		assert(bbstop==h.bbstop) : s+"\n"+s2;
//		assert(bbchrom==h.bbchrom) : s+"\n"+s2;
//		
//		assert(rname.equals(h.rname)) : s+"\n"+s2;
	}
	
	public CustomHeader(Read r){
		
		rnum=r.pairnum();
		id=r.numericID;//(rnum==1 && r.mate!=null ? r.mate.numericID : r.numericID);
//		assert(rnum==0) : id;
		bbchrom=r.chrom;
		bbstart=r.start;
		int bbstop=r.stop;
		match=r.match;
		strand=r.strand();

		start=0;
		stop=0;
		rname=null;
		int reflen=2000000000;
		
		if(Data.GENOME_BUILD>=0){
			int bbscaffold=Data.scaffoldIndex(bbchrom, (bbstart+bbstop)/2);
			byte[] name1=Data.scaffoldNames[bbchrom][bbscaffold];
			start=Data.scaffoldRelativeLoc(bbchrom, bbstart, bbscaffold);
			stop=start-bbstart+bbstop;
			rname=new String(name1);
			reflen=Data.scaffoldLengths[bbchrom][bbscaffold];
		}
		insert=(r.mate==null ? 0 : r.insertSizeMapped(false));
		

		
//		String s=toString();
//		Header h=new Header("SYN_"+s, rnum);
//		String s2=h.toString();
//		assert(rnum==h.rnum) : s+"\n"+s2;
//		assert(id==h.id) : s+"\n"+s2;
//		assert(start==h.start) : s+"\n"+s2;
//		assert(stop==h.stop) : s+"\n"+s2;
//		assert(insert==h.insert) : s+"\n"+s2;
//		assert(strand==h.strand) : s+"\n"+s2;
//		
//		assert(bbstart==h.bbstart) : s+"\n"+s2;
//		assert(bbstop==h.bbstop) : s+"\n"+s2;
//		assert(bbchrom==h.bbchrom) : s+"\n"+s2;
//		
//		assert(rname.equals(h.rname)) : rname+"\n"+h.rname;
	}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder(64);
		appendTo(bb);
		return bb.toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){
//		if(rnum==0){bb.append("SYN_");}
//		bb.append(rnum);
		bb.append(id);
		bb.append('_').append(start);
		bb.append('_').append(stop);
		bb.append('_').append(insert);
		bb.append('_').append(Shared.strandCodes2[strand]);
		bb.append('_').append(bbstart);
		bb.append('_').append(bbchrom);
		bb.append('_').append(match==null ? "." : new String(match));
		bb.append('_').append(rname==null ? "." : encodeRname(rname));
		return bb;
	}
	
	public static String customID_old(Read r){
		if(!FASTQ.TAG_CUSTOM){return r.id;}
		
		if(FASTQ.DELETE_OLD_NAME){
			assert(false) : "Seems odd so I added this assertion.  I don't see anywhere it was being used. Use -da flag to override.";
			r.id=null;
		}

		ByteBuilder sb=new ByteBuilder();
		if(r.id==null /*|| DELETE_OLD_NAME*/){
			sb.append(r.numericID);
		}else{
			sb.append(r.id);
		}
		if(r.chrom>-1 && r.stop>-1){
			if(FASTQ.TAG_CUSTOM_SIMPLE){
				sb.append('_');
				sb.append(r.strand()==0 ? '+' : '-');
			}else{
				sb.append("_chr");
				sb.append(r.chrom);
				sb.append('_');
				sb.append((int)r.strand());
				sb.append('_');
				sb.append(r.start);
				sb.append('_');
				sb.append(r.stop);
			}

			if(Data.GENOME_BUILD>=0){
				final int chrom1=r.chrom;
				final int start1=r.start;
				final int stop1=r.stop;
				final int idx1=Data.scaffoldIndex(chrom1, (start1+stop1)/2);
				final byte[] name1=Data.scaffoldNames[chrom1][idx1];
				final int a1=Data.scaffoldRelativeLoc(chrom1, start1, idx1);
				final int b1=a1-start1+stop1;
				sb.append('_');
				sb.append(a1);
				if(FASTQ.TAG_CUSTOM_SIMPLE){
					sb.append('_');
					sb.append(b1);
				}
				sb.append('_');
				sb.append(new String(name1));
			}
		}
		
		if(FASTQ.ADD_PAIRNUM_TO_CUSTOM_ID){
			sb.append(' ');
			sb.append(r.pairnum()+1);
			sb.append(':');
		}else if(FASTQ.ADD_SLASH_PAIRNUM_TO_CUSTOM_ID){
			if(FASTQ.SPACE_SLASH){sb.append(' ');}
			sb.append('/');
			sb.append(r.pairnum()+1);
		}
		
		return sb.toString();
	}
	
	public static String toString(Read r){
		CustomHeader h1=new CustomHeader(r);
		CustomHeader h2=(r.mate==null ? null : new CustomHeader(r.mate));
		final String s;
		if(r.pairnum()==0){
			s=toString(h1, h2, r.pairnum());
		}else{
			s=toString(h2, h1, r.pairnum());
		}
//		System.err.println(r.pairnum());
		return s;
	}
	
	public static String toString(CustomHeader h1, CustomHeader h2, int rnum){
		ByteBuilder bb=new ByteBuilder();
		bb.append("SYN_");
		h1.appendTo(bb);
		if(h2!=null){
			bb.append('&');
			h2.appendTo(bb);
		}
		
		if(FASTQ.ADD_PAIRNUM_TO_CUSTOM_ID){
			bb.append(' ');
			bb.append(rnum+1);
			bb.append(':');
		}else if(FASTQ.ADD_SLASH_PAIRNUM_TO_CUSTOM_ID){
			if(FASTQ.SPACE_SLASH){bb.append(' ');}
			bb.append('/');
			bb.append(rnum+1);
		}
		return bb.toString();
	}
	
	public static String decodeRname(String rname){
		ByteBuilder bb=new ByteBuilder(rname.length());
		char prev='.';
		for(int i=0; i<rname.length(); i++){
			char c=rname.charAt(i);
			if(prev=='!'){
				bb.append(escapeArray[c]);
				prev='.';
			}else if(c!='!'){
				bb.append(decodeArray[c]);
				prev=c;
			}else{
				prev=c;
			}
		}
		return bb.toString();
	}
	
	public static String encodeRname(String rname){
		
		int found=0;
		for(int i=0; i<rname.length(); i++){
			char c=rname.charAt(i);
			if(c=='!' || c=='&' || c==' ' || c=='$' || c=='{'){found++;}
		}
		if(found<1){return rname;}
		
		ByteBuilder bb=new ByteBuilder(rname.length()+found);
		for(int i=0; i<rname.length(); i++){
			char c=rname.charAt(i);
			if(c=='!'){bb.append("!!");}
			else if(c==' '){bb.append("$");}
			else if(c=='_'){bb.append("{");}
			else if(c=='&'){bb.append("!a");}
			else if(c=='$'){bb.append("!d");}
			else if(c=='{'){bb.append("!b");}
			else{bb.append(c);}
		}
		return bb.toString();
	}
	
	public static int getPairnum(String s){
		int idx=s.lastIndexOf(' ');
		if(idx<0){idx=s.lastIndexOf('/');}
		assert(idx>0) : s;
		idx++;
		if(s.charAt(idx)=='/'){idx++;}
//		assert(false) : s.substring(idx-1);
		return s.charAt(idx)-'1';
	}

	public int rnum;
	public long id;
	public int start;
	public int stop;
	public int insert;
	public int strand;
	public int bbstart;
	public int bbstop(){return bbstart-start+stop;}
	public int bbchrom;
	public byte[] match;
	public String rname;

	private static final char ESCAPE='!';
	private static final char MIDDLE='&';

	private static final byte[] decodeArray;
	private static final byte[] escapeArray;
	
	static{
		decodeArray=new byte[128];
		escapeArray=new byte[128];
		for(int i=0; i<128; i++){
			decodeArray[i]=escapeArray[i]=(byte)i;
			if(i=='$'){decodeArray[i]=' ';}
			else if(i=='{'){decodeArray[i]='_';}
			else if(i=='a'){escapeArray[i]='&';}
			else if(i=='b'){escapeArray[i]='{';}
			else if(i=='d'){escapeArray[i]='$';}
		}
	}
	
}
