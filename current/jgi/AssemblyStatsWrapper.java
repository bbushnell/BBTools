package jgi;

import java.io.File;
import java.util.ArrayList;

import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Apr 17, 2013
 *
 */
public class AssemblyStatsWrapper {
	
	public static void main(String[] args){
		ArrayList<String> alist=new ArrayList<String>();
		ArrayList<String> ilist=new ArrayList<String>();
		
		alist.add("");
		alist.add("header=t");
		alist.add("showspeed=f");
		alist.add("addname=t");
		alist.add("k=0");
		
		for(String arg : args){
			if(!arg.contains("=") && Tools.canRead(arg)){
				ilist.add("in="+arg);
			}else{
				String[] split=arg.split("=");
				if(split[0].equalsIgnoreCase("in") || split[0].equalsIgnoreCase("ref")){
					if(split.length>1){
						if(new File(split[1]).exists()){
							ilist.add(arg);
						}else{
							String[] split2=split[1].split(",");
							for(String s : split2){
								ilist.add("in="+s);
							}
						}
					}
				}else{
					alist.add(arg);
				}
			}
		}
		
		String[] args2=alist.toArray(new String[alist.size()]);
		for(int i=0; i<ilist.size(); i++){
			String s=ilist.get(i);
//			System.err.println("Processing "+s);
			args2[0]=s;
			if(i>0){
				args2[1]="header=f";
//				AssemblyStats2.reset();
				System.gc();
				synchronized(AssemblyStatsWrapper.class){
					try {
						AssemblyStatsWrapper.class.wait(100);
					} catch (InterruptedException e) {}
				}
				Thread.yield();
			}
			AssemblyStats2 as2=new AssemblyStats2(args2);
			if(i>0){
				AssemblyStats2.overwrite=false;
				AssemblyStats2.append=true;
			}
			as2.process();
		}
		
	}
	
}
