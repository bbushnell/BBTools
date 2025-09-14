package gff;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

public class GbffFeature {

	public GbffFeature(final ArrayList<byte[]> lines0, final String typeString, final String accessionString){
		accession=accessionString;
		setType(typeString);
		parseSlow(lines0);
		if(type==rRNA){
			setSubtype();
		}
		if(stop<start){error=true;}
	}
	
	private void parseSlow(final ArrayList<byte[]> lines0){
		ArrayList<byte[]> lines=fixLines(lines0);
		parseStartStop(lines.get(0));
		for(int i=1; i<lines.size(); i++){
			byte[] line=lines.get(i);
			if(Tools.startsWith(line, "product=")){
				product=parseLine(line);
			}else if(Tools.startsWith(line, "locus_tag=")){
				locus_tag=parseLine(line);
			}else if(Tools.equals(line, "pseudo")){
				pseudo=true;
			}
			
//			else if(Tools.startsWith(line, "ID=")){
//				id=parseLine(line);
//			}else if(Tools.startsWith(line, "Name=")){
//				name=parseLine(line);
//			}
		}
//		System.err.println("\nvvvvv");
//		for(byte[] line : lines0){
//			System.err.println("'"+new String(line)+"'");
//		}
//		for(byte[] line : lines){
//			System.err.println("'"+new String(line)+"'");
//		}
//		System.err.println("^^^^^");
	}
	
	ArrayList<byte[]> fixLines(ArrayList<byte[]> lines){
		ArrayList<byte[]> fixed=new ArrayList<byte[]>();
		ByteBuilder bb=new ByteBuilder();
		for(byte[] line : lines){
			if(bb.length()>0 && line[21]=='/'){
				fixed.add(bb.toBytes());
				bb.clear();
			}
			append(bb, line);
		}
		if(bb.length()>0){
			fixed.add(bb.toBytes());
			bb.clear();
		}
		return fixed;
	}
	
	void append(ByteBuilder bb, byte[] line){
		assert(line[20]==' ');
		assert(line.length>21);
//		assert(line[21]!=' ') : "'"+new String(line)+"'";
		if(line[21]=='/'){
			bb.append(line, 22, line.length-22);
		}else{
//			System.err.println(line.length+", "+21+", "+(line.length-21+1)+"\n'"+new String(line)+"'");
			if(bb.length>0){bb.append(' ');}
			bb.append(line, 21, line.length-21);
		}
	}
	
	void setType(String typeString){
		int x=Tools.find(typeString, typeStrings);
		assert(x>=0) : x+", "+typeString;
		type=x;
	}
	
	void parseStartStop(final byte[] line0){
		byte[] line=line0;
		
		if(line[0]=='c'){
			assert(Tools.startsWith(line, "complement("));
			line=Arrays.copyOfRange(line, 11, line.length-1);
			strand=Shared.MINUS;
		}
		if(line[0]=='j'){
			assert(Tools.startsWith(line, "join("));
			line=Arrays.copyOfRange(line, 5, line.length-1);
			strand=Shared.MINUS;
		}
		
		int i=0;
		for(start=0; i<line.length; i++){
			int x=line[i];
			if(x=='.'){break;}
			else if(x!='<'){
				if(Tools.isDigit(x)){
					start=start*10+(x-'0');
				}else{
					//if(!error){System.err.println(new String(line0)+"\n"+new String(line));}
					error=true;
				}
			}
		}
//		while(line[i]=='.'){i++;} //Not needed
		for(stop=0; i<line.length; i++){
			int x=line[i];
			if(x=='.' || x==','){
				stop=0;
			}else if(x==' '){
				//do nothing; line wrap
			}else if(x!='>'){
				if(Tools.isDigit(x)){
					stop=stop*10+(x-'0');
				}else{
					//if(!error){System.err.println(new String(line0)+"\n"+new String(line));}
					error=true;
				}
			}
		}
	}
	
	String parseLine(byte[] line){
		String[] split=Tools.equalsPattern.split(new String(line));
		String s=split[1];
		return s.substring(1, s.length()-1);
	}
	
	void setSubtype(){
		subtype=-1;
		if(product==null){return;}
		String[] split=Tools.spacePattern.split(product);
		subtype=Tools.find(split[0], typeStrings);
//		assert(false) : type+", "+subtype+", "+split[0]+", "+this.toString()+"\n"+product;
	}
	
	public void toGff(ByteStreamWriter bsw) {
		ByteBuilder bb=bsw.getBuffer();
		appendGff(bb);
		bb.nl();
		bsw.flushBuffer(false);
	}
	
	public ByteBuilder appendGff(ByteBuilder bb) {
//		bsw.print("#seqid	source	type	start	end	score	strand	phase	attributes\n".getBytes());
		bb.append(accession).tab();
		bb.append('.').tab();
		bb.append((pseudo && type==GENE) ? "pseudogene" : typeStringsGff[type]).tab();
		bb.append(start).tab();
		bb.append(stop).tab();
		bb.append('.').tab();
		bb.append(Shared.strandCodes2[strand]).tab();
		bb.append('.').tab();
		
		boolean attributes=false;
//		if(id!=null){
//			bb.append("ID=").append(id);
//			attributes=true;
//		}
//		if(name!=null){
//			if(attributes){bb.append(';');}
//			bb.append("Name=").append(name);
//			attributes=true;
//		}
		if(product!=null){
			if(attributes){bb.append(';');}
			bb.append("product=").append(product);
			attributes=true;
		}
		if(locus_tag!=null){
			if(attributes){bb.append(';');}
			bb.append("locus_tag=").append(locus_tag);
			attributes=true;
		}
		if(subtype>-1){
			if(attributes){bb.append(';');}
			bb.append("subtype=").append(typeStringsGff[subtype]);
			attributes=true;
		}
		if(!attributes){bb.append('.');}
		return bb;
	}
	
	
	@Override
	public String toString(){
		return appendGff(new ByteBuilder()).toString();
	}

	public int type=-1;
	public int subtype=-1;
	//TODO: could have coding amino, for tRNA
	public String product;
	public String locus_tag;
//	public String id;
//	public String name;
	
	public int start;
	public int stop;
	public byte strand=Shared.PLUS;
	public String accession;
	public boolean pseudo=false;
	public boolean error=false;

	public static final String[] typeStrings={"gene", "CDS", "rRNA", "tRNA", "ncRNA", "repeat_region", 
			"5'UTR", "3'UTR", "intron", "exon", "5S", "16S", "23S"};
	public static final String[] typeStringsGff={"gene", "CDS", "rRNA", "tRNA", "ncRNA", "repeat_region", 
			"five_prime_UTR", "three_prime_UTR", "intron", "exon", "5S", "16S", "23S"};
	
	//types
	public static final int GENE=0, CDS=1, rRNA=2, tRNA=3, ncRNA=4, repeat_region=5, UTR5=6, UTR3=7, intron=8, exon=9;
	//subtypes
	public static final int r5S=10, r16S=11, r23S=12;
	
}
