package tax;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Timer;
import shared.Tools;

public class ImgRecord implements Serializable {
	
	private static final long serialVersionUID = 6438551103300423985L;
	
	public static void main(String[] args){
		String in=args[0];
		String out=args.length>1 ? args[1] : null;

		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testOutputFiles(true, false, false, out)){
			throw new RuntimeException("\nCan't write to some output files.\n");
		}
		Timer t=new Timer();
		HashMap<Long, ImgRecord> map=toMap(in, TaxTree.IMG_HQ);
		t.stop();
		System.err.println(map.size()+"; "+t);
//		if(out!=null){ReadWrite.writeObjectInThread(map, out, false);}
		if(out!=null){writeAsText(map, out);}
	}
	
	private static void writeAsText(HashMap<Long, ImgRecord> map, String out){
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		for(Entry<Long, ImgRecord> e : map.entrySet()){
			tsw.println(e.toString());
		}
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append(imgID);
		sb.append('\t').append(taxID);
		sb.append('\t').append(name);
		return sb.toString();
	}
	
	public static HashMap<Long, ImgRecord> toMap(String fname, boolean highQuality){
		ImgRecord[] array=toArray(fname, highQuality);
		HashMap<Long, ImgRecord> map=new HashMap<Long, ImgRecord>((3+array.length*4)/3);
		for(ImgRecord ir : array){
			map.put(ir.imgID, ir);
		}
		return map;
	}
	
	public static ImgRecord[] toArray(String fname, boolean highQuality){
		TextFile tf=new TextFile(fname, false);
		ArrayList<ImgRecord> list=new ArrayList<ImgRecord>();
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.length()<1 || !Tools.isDigit(line.charAt(0))){
				//do nothing
			}else{
				ImgRecord record=new ImgRecord(line);
				if(!highQuality || record.highQuality){list.add(record);}
			}
		}
		tf.close();
		return list.toArray(new ImgRecord[0]);
	}
	
	public ImgRecord(String line){
		String[] split=line.split("\t");
		
		imgID=Long.parseLong(split[0]);
		name=(storeName ? split[1] : null);
		try {
			taxID=(split[2]==null || split[2].length()<1 ? -1 : Integer.parseInt(split[2]));
		} catch (NumberFormatException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			System.err.println(line);
			throw new RuntimeException();
		}
		isPublic=Parse.parseYesNo(split[3]);
		obsolete=Parse.parseYesNo(split[4]);
		genomeType=find(split[5], typeArray);
		boolean hq=false;
		if(split.length>7){
			try {
				hq=Parse.parseYesNo(split[7]);
			} catch (Exception e) {
				System.err.println(Arrays.toString(split));
				assert(false);
			}
		}
		highQuality=hq;
	}
	
	public final long imgID;
	public final int taxID;
	public final int genomeType;
	public final boolean isPublic;
	public final boolean obsolete;
	public final boolean highQuality;
	public final String name;
	public final String path(){return "/global/dna/projectdirs/microbial/img_web_data/taxon.fna/"+imgID+".fna";}
	
	final int ISOLATE=0, SINGLE_CELL=1, METAGENOME=2;
	final String[] typeArray={"isolate", "single_cell", "metagenome"};
	private static int find(String s, String[] array){
		for(int i=0; i<array.length; i++){
			if(array[i].equals(s)){return i;}
		}
		return -1;
	}
	
	public static boolean storeName=true;
	public static HashMap<Long, ImgRecord> imgMap;
//	public static final String DefaultDumpFile="/global/cfs/cdirs/bbtools/tax/imgTaxDump.txt.gz";
	public static final String DefaultDumpFile="/global/u1/i/img/adhocDumps/taxonDumpForBrian.txt";
	
}
