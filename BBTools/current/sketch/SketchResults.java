package sketch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import fileIO.TextStreamWriter;
import json.JsonObject;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntHashMap;
import tax.TaxTree;

public class SketchResults extends SketchObject {
	
	SketchResults(Sketch s){
		sketch=s;
	}
	
	SketchResults(Sketch s, ArrayList<Sketch> sketchList_, int[][] taxHits_){
		sketch=s;
		refSketchList=sketchList_;
		taxHits=taxHits_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void addMap(ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, CompareBuffer buffer) {

		if(map.isEmpty()){return;}
		list=addToList(map, params, list);
		
		if((true || params.needContamCounts())){
			recompare(buffer, params);
		}
	}
	
	public void recompare(CompareBuffer buffer, DisplayParams params){
//		assert(makeIndex || !AUTOSIZE);
		
		assert(!sketch.merged());
		sketch.mergeBitSets();
		
//		System.err.println(sketch.compareBitSet());
//		assert(false) : sketch.compareBitSet().getClass();
		
		for(Comparison c : list){
			c.recompare(buffer, taxHits, params.contamLevel());
		}
		Collections.sort(list, params.comparator);
		Collections.reverse(list);
	}
	
	private static ArrayList<Comparison> addToList(ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, ArrayList<Comparison> old){
		
//		System.err.println(map.size());
//		System.err.println(map.keySet());

//		final TaxFilter white=params.taxFilterWhite;
//		final TaxFilter black=params.taxFilterBlack;
//		final boolean noFilter=(white==null && black==null);
		final int size=map.size();
		ArrayList<Comparison> al=(old==null ? new ArrayList<Comparison>(size) : old);
		for(Entry<Integer, Comparison> e : map.entrySet()){
			final Comparison c=e.getValue();
			al.add(c);
//			if(noFilter || c.passesFilter(white, black)){
//				al.add(c);
//			}
		}
		Shared.sort(al, params.comparator);
		Collections.reverse(al);
		
		//Apply records per level filter
		if(params.recordsPerLevel>0 && al.size()>params.recordsPerLevel && al.get(0).hasQueryTaxID()){
			int[] count=new int[TaxTree.numTaxLevelNamesExtended];
			int removed=0;
			for(int i=0; i<al.size(); i++){
				Comparison c=al.get(i);
				int calevel=c.commonAncestorLevelInt();
				count[calevel]++;
				if(count[calevel]>params.recordsPerLevel){
					al.set(i, null);
					removed++;
				}
			}
			if(removed>0){Tools.condenseStrict(al);}
		}
		
		final long limit=(params.maxRecords*2+5);
		while(al.size()>limit){
			al.remove(al.size()-1);
		}
		return al;
	}
	
	public boolean isEmpty(){
		return list==null || list.isEmpty();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Tax Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public int primaryTax(int level){
		//I have no idea how to implement this...
		IntHashMap map=new IntHashMap();
		assert(false);
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Print Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static String recordBreak="\n"; //"\n\n"
	
	void writeResults(DisplayParams params, TextStreamWriter tsw){
		ByteBuilder sb=toText(params);
		tsw.print(sb);
	}
	
	public ByteBuilder toText(DisplayParams params){
		assert(params.postParsed);
		if(sketch.hasSSU()){
			if(params.comparator==Comparison.SSUComparator){
				alignSSUs(params.maxRecords*4);//This should be enough...
				list.sort(params.comparator);
				Collections.reverse(list);
			}else if(params.printSSU()){
				alignSSUs(params.maxRecords);
			}
		}
		if(params.json()){
			JsonObject j=params.toJson(this);
			return j.toText();
		}
		final ByteBuilder sb=params.queryHeader(sketch);
		if(params.format==DisplayParams.FORMAT_QUERY_REF_ANI || params.format==DisplayParams.FORMAT_CONSTELLATION){
			if(list==null || list.isEmpty()){return sb;}
			int idx=0;
			int prevTaxID=0;
			for(Comparison c : list){
				assert(!params.printSSU() || !c.needsAlignment());
				params.formatComparison(c, sb, prevTaxID);
				prevTaxID=c.taxID();
				idx++;
				if(idx>=params.maxRecords){break;}
			}
		}else{
			sb.append(recordBreak);

			if(list==null || list.isEmpty()){
				sb.append("No hits.\n");
			}else{
				if(params.format==DisplayParams.FORMAT_MULTICOLUMN){sb.append(params.header()).append('\n');}
				int idx=0;
				int prevTaxID=0;
				for(Comparison c : list){
					params.formatComparison(c, sb, prevTaxID);
					prevTaxID=c.taxID();
					idx++;
					if(idx>=params.maxRecords){break;}
				}
			}
		}
		return sb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Alignment           ----------------*/
	/*--------------------------------------------------------------*/
	
	void alignSSUs(int maxRecords){
		if(!sketch.hasSSU()){return;}
//		if(list!=null && list.size()>0){
//			int idx=0;
//			for(Comparison c : list){
//				c.ssuIdentity();
//				idx++;
//				if(idx>=maxRecords){break;}
//			}
//		}
		assert(alignerPool!=null);
		alignerPool.addJobs(list, maxRecords);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final Sketch sketch;
	public ArrayList<Sketch> refSketchList;
	public int[][] taxHits;
	public ArrayList<Comparison> list;
	public int totalRecords=0;
	
}
