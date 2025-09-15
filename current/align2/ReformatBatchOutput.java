package align2;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.TextFile;
import shared.Tools;

public class ReformatBatchOutput {
	
//	Elapsed:	31.7
//
//	Mapping Statistics for 0s_default.sam:
//	mapped:                	100.00%
//	retained:              	96.06%
//	discarded:             	0.00%
//	ambiguous:             	3.94%
//
//	Strict correctness (both ends exactly correct):
//	true positive:         	96.06%
//	false positive:        	0.00%
//
//	Loose correctness (one end approximately correct):
//	true positive:         	96.06%
//	false positive:        	0.00%
//
//	false negative:        	0.00%
//	Elapsed:	2.34
//	Elapsed:	20.51
	
	
//	Elapsed:	0.33
//
//	Mapping Statistics for bwa_0S_0I_0D_0U_0N_r100.sam:
//	primary alignments:    	100 found of 100 expected
//	secondary alignments:  	0 found
//	mapped:                	100.000%
//	retained:              	97.000%
//	discarded:             	 0.000%
//	ambiguous:             	 3.000%
//
//	Strict correctness (both ends exactly correct):
//	true positive:         	97.000%
//	false positive:        	 0.000%
//
//	Loose correctness (one end approximately correct):
//	true positive:         	97.000%
//	false positive:        	 0.000%
//
//	false negative:        	 0.000%
	
	
	public static void main(String[] args){
		TextFile tf=new TextFile(args[0], false);
		String[] lines=tf.toStringLines();
		ArrayList<String> list=new ArrayList<String>();
		
		int mode=0;
		
		System.out.println(header());
		
		for(String s : lines){
			if(s.startsWith("Elapsed:")){
				if(!list.isEmpty()){
					process(list); //failure
					list.clear();
					mode=0;
				}
				mode++;
			}
			
			if(mode>0){
				list.add(s);
				if(s.startsWith("false negative:")){
					process(list);
					list.clear();
					mode=0;
				}
			}
		}
	}
	
	
	public static String header() {
		return("program\tfile\tvartype\tcount\treads\tprimary\tsecondary\ttime\tmapped\tretained\tdiscarded\tambiguous\ttruePositive\t" +
				"falsePositive\ttruePositiveL\tfalsePositiveL\tfalseNegative");
	}
	
	//bwa_1S_0I_0D_0U_0N_r400000x100.sam
	public static int getReads(String name){
//		String[] split=name.substring(0, name.length()-4).split("_");
		String[] split=name.split("_");
		String r=(split[split.length-1]);
		if(r.charAt(0)=='r' && Tools.isDigit(r.charAt(r.length()-1))){
			assert(r.charAt(0)=='r') : Arrays.toString(split)+", "+name;
			r=r.substring(1);
			if(r.contains("x")){
				r=r.substring(0, r.indexOf('x'));
			}
			return Integer.parseInt(r);
		}else{
			for(String s : split){
				if(s.endsWith("bp") && s.contains("x") && Tools.isDigit(s.charAt(0))){
					r=s.substring(0, s.indexOf('x')-1);
					return Integer.parseInt(r);
				}
			}
		}
		return 0;
	}
	
	public static char getVarType(String name){
//		String[] split=name.substring(0, name.length()-4).split("_");
		String[] split=name.split("_");
		for(String s : split){
			char c=s.charAt(0);
			if(Tools.isDigit(c) && c!='0' && !s.endsWith("bp")){
				return s.charAt(s.length()-1);
			}
		}
		return '?';
	}
	
	public static int getCount(String name){
//		String[] split=name.substring(0, name.length()-4).split("_");
		String[] split=name.split("_");
		for(String s : split){
			char c=s.charAt(0);
			if(Tools.isDigit(c) && c!='0' && !s.endsWith("bp")){
				String r=s.substring(0, s.length()-1);
				return Integer.parseInt(r);
			}
		}
		return 0;
	}
	
	public static String getProgram(String name){
		return name.substring(0, name.indexOf('_'));
	}
	
	

	public static void process(ArrayList<String> list){
		
		String name=null;
//		String count=null;
		String time=null;
		StringBuilder sb=new StringBuilder();
		
		int primary=0;
		int secondary=0;
		int expected=0;
		
		for(String s : list){
			String[] split=s.split("\t");
			String a=split[0];
			String b=split.length>1 ? split[1] : null;
			if(a.equals("Elapsed:")){
				time=b;
			}else if(a.startsWith("lines:")){
				//do nothing
			}else if(a.startsWith("Mapping Statistics for ")){
				name=a.replace("Mapping Statistics for ", "").replace(".sam:", "");
			}else if(a.startsWith("primary alignments:")){
				assert(b!=null) : "Bad line: "+s;
				b=b.replace(" found of ", "_");
				b=b.replace(" expected", "");
				String[] split2=b.split("_");
				primary=Integer.parseInt(split2[0]);
				expected=Integer.parseInt(split2[1]);
			}else if(a.startsWith("secondary alignments:")){
				assert(b!=null) : "Bad line: "+s;
				b=b.replace(" found", "");
				secondary=Integer.parseInt(b);
			}else if(b!=null){
				assert(!b.contains("found")) : "\na='"+a+"'\nb='"+b+"'\n"+a.equals("primary alignments:");
				sb.append('\t').append(b.replace("%", ""));
			}
			
		}
		
//		if(name!=null){
//			count="";
//
//			String[] split=name.split("_");
//
//			for(String s : split){
//				if(s!=null && s.length()>0 && s.charAt(0)!='0'){
//					for(int i=0; i<s.length(); i++){
//						char c=s.charAt(i);
//						if(Tools.isDigit(c)){count=count+c;}
//						else{break;}
//					}
//				}
//				if(count.length()>0){break;}
//			}
//		}
		
		String prg=null;
		char type='S';
		int reads=1;
		int vars=0;
		
		if(name!=null){
			try {
				prg=getProgram(name);
				type=getVarType(name);
				reads=getReads(name);
				vars=getCount(name);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		System.out.println(prg+"\t"+name+"\t"+type+"\t"+vars+"\t"+reads+"\t"+primary+"\t"+secondary+"\t"+time+sb);
		
	}
	
}
