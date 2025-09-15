package gff;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;

public class GbffFile {
	
	public static void main(String[] args){
		String gbff=args[0];
		String gff=(args.length>1 ? args[1] : "stdout.gff");

		if(gbff.indexOf('=')>=0){gbff=gbff.split("=")[1];}
		if(gff.indexOf('=')>=0){gff=gff.split("=")[1];}
		
		FileFormat ffin=FileFormat.testInput(gbff, ".gbff", true);
		FileFormat ffout=FileFormat.testOutput(gff, FileFormat.GFF, null, true, true, false, false);
		GbffFile file=new GbffFile(ffin);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		file.toGff(bsw, true);
		bsw.poisonAndWait();
	}
	
	
//	##gff-version 3
//	#!gff-spec-version 1.21
//	#!processor NCBI annotwriter
//	#!genome-build IMG-taxon 2724679794 annotated assembly
//	#!genome-build-accession NCBI_Assembly:GCF_900182635.1
//	#!annotation-date 07/14/2019 01:52:19
//	#!annotation-source NCBI RefSeq 
//	##sequence-region NZ_FXTD01000001.1 1 528269
//	##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=413815
	
	public void toGff(ByteStreamWriter bsw, boolean printHeader){
		if(printHeader){
			bsw.println("##gff-version 3".getBytes());
			bsw.println(("#BBTools "+Shared.BBTOOLS_VERSION_STRING+" GbffToGff").getBytes());
			bsw.println("#seqid	source	type	start	end	score	strand	phase	attributes".getBytes());
		}
		for(GbffLocus locus=nextLocus(); locus!=null; locus=nextLocus()){
			locus.toGff(bsw);
		}
	}
	
	public GbffFile(FileFormat ff_) {
		ff=ff_;
		assert(ff.format()==FileFormat.GBFF) : ff;
		reset();
	}
	
	public synchronized void reset(){
		if(bf!=null){
			bf.close();
			bf=null;
		}
		bf=ByteFile.makeByteFile(ff, FileFormat.GBFF);
		line=bf.nextLine();
		if(line==null){bf.close();}//empty
	}
	
	public GbffLocus nextLocus(){
		assert(bf!=null);
		if(line==null){return null;}
		assert(Tools.startsWith(line, "LOCUS ")) : "Expecting: 'LOCUS ...'\nGot: '"+new String(line)+"'";
		ArrayList<byte[]> lines=new ArrayList<byte[]>();
		lines.add(line);
		boolean sequence=false;
		for(line=bf.nextLine(); line!=null && (line.length==0 || line[0]!='L' || !Tools.startsWith(line, "LOCUS ")); line=bf.nextLine()){
			if(line.length>0){
				final byte b=line[0];
				if(b=='/'){
					//skip
				}else if(b=='O' && Tools.startsWith(line, "ORIGIN ")){
					sequence=true;
				}else if(b==' ' && sequence){
						//do nothing
				}else{
					sequence=false;
					lines.add(line);
				}
			}
		}
		if(line==null){bf.close();}
		return new GbffLocus(lines);
	}
	
	private final FileFormat ff;
	private ByteFile bf;
	private byte[] line=null;
	
}
