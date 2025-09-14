package tax;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Timer;
import shared.Tools;

public class ImgRecord2 implements Serializable {
	
	private static final long serialVersionUID = 8596055205924485293L;
	
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
		HashMap<Long, ImgRecord2> map=toMap(in);
		t.stop();
		System.err.println(map.size()+"; "+t);
		if(out!=null){ReadWrite.writeObjectInThread(map, out, false);}
	}
	
	public static HashMap<Long, ImgRecord2> toMap(String fname){
		ImgRecord2[] array=toArray(fname);
		HashMap<Long, ImgRecord2> map=new HashMap<Long, ImgRecord2>((3+array.length*4)/3);
		for(ImgRecord2 ir : array){
			map.put(ir.imgID, ir);
		}
		return map;
	}
	
	public static ImgRecord2[] toArray(String fname){
		TextFile tf=new TextFile(fname, false);
		ArrayList<ImgRecord2> list=new ArrayList<ImgRecord2>();
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.length()<1 || !Tools.isDigit(line.charAt(0))){
				//do nothing
			}else{
				ImgRecord2 record=new ImgRecord2(line);
				list.add(record);
			}
		}
		tf.close();
		return list.toArray(new ImgRecord2[0]);
	}
	
	public ImgRecord2(String line){
		String[] split=line.split("\t");
		
		imgID=Long.parseLong(split[0]);
		taxID=(split[1]==null || split[1].length()<1 ? 0 : Integer.parseInt(split[1]));
		name=(storeName ? split[2] : null);
	}
	
	public static final long parseImgId(String header, boolean doAssertions){
		final char delimiter=TaxTree.ncbiHeaderDelimiter(header);
		assert(!doAssertions || delimiter!=' ') : header;
//		System.err.println("A: '"+delimiter+"'");
		if(delimiter==' '){return -1;}
		
		final int idx=header.indexOf("img"+delimiter);
//		System.err.println("B: "+idx);
		assert(!doAssertions || idx>=0) : "Could not img id number from "+header;
		if(idx<0){return -1;}
		
		long img=0;
		for(int i=idx+4; i<header.length(); i++){
			final char c=header.charAt(i);
			if(c==delimiter || c==' '){break;}
			assert(Tools.isDigit(c)) : c+", '"+header+"'";
			img=img*10+(c-'0');
		}
		
		assert(!doAssertions || img>0) : "Could not img id number from "+header;
		return img>0 ? img : -1;
	}
	
	public static boolean storeName=true;
	public final long imgID;
	public final int taxID;
	public final String name;
	
}
