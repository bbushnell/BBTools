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

/**
 * Manages and processes results from genomic sketch comparisons.
 * Provides comprehensive result tracking, filtering, and analysis capabilities
 * for sketch-based taxonomic classification and sequence similarity evaluation.
 * Supports concurrent result collection, dynamic result list management,
 * and multi-level taxonomic hit tracking.
 *
 * @author Brian Bushnell
 */
public class SketchResults extends SketchObject {
	
	/** Constructs a SketchResults object with a single sketch.
	 * @param s The query sketch to associate with this result set */
	SketchResults(Sketch s){
		sketch=s;
	}
	
	/**
	 * Constructs a SketchResults object with sketch and reference data.
	 * @param s The query sketch
	 * @param sketchList_ List of reference sketches for comparison
	 * @param taxHits_ Two-dimensional array tracking taxonomic hits
	 */
	SketchResults(Sketch s, ArrayList<Sketch> sketchList_, int[][] taxHits_){
		sketch=s;
		refSketchList=sketchList_;
		taxHits=taxHits_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds comparison results from a concurrent map to the result list.
	 * Converts the concurrent hashmap of comparisons into a sorted list
	 * and triggers recomparison if contamination counts are needed.
	 *
	 * @param map Concurrent map containing sketch comparisons keyed by ID
	 * @param params Display parameters controlling comparison behavior
	 * @param buffer Comparison buffer for recomparison calculations
	 */
	public void addMap(ConcurrentHashMap<Integer, Comparison> map, DisplayParams params, CompareBuffer buffer) {

		if(map.isEmpty()){return;}
		list=addToList(map, params, list);
		
		if((true || params.needContamCounts())){
			recompare(buffer, params);
		}
	}
	
	/**
	 * Recompares all results using merged bit sets for more accurate metrics.
	 * Merges the query sketch's bit sets, then recomputes comparison metrics
	 * for all comparisons in the result list. Sorts results by the specified
	 * comparator after recomparison.
	 *
	 * @param buffer Comparison buffer providing computational resources
	 * @param params Display parameters containing contamination level settings
	 */
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
	
	/**
	 * Converts concurrent map to sorted comparison list with level filtering.
	 * Applies records-per-level filtering to limit results at each taxonomic level
	 * and truncates the list to maximum record limits. Maintains sorting based
	 * on the provided comparator.
	 *
	 * @param map Concurrent map of comparisons to convert
	 * @param params Display parameters controlling filtering and limits
	 * @param old Existing comparison list to append to (may be null)
	 * @return Filtered and sorted list of comparisons
	 */
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
	
	/** Checks if the result list is empty or null.
	 * @return true if no comparison results are present, false otherwise */
	public boolean isEmpty(){
		return list==null || list.isEmpty();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Tax Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * FunctionUnclear: Intended to determine primary taxonomic ID at given level.
	 * Current implementation contains assertion failure indicating incomplete
	 * development. Method signature suggests taxonomic level processing but
	 * implementation is not functional.
	 *
	 * @param level Taxonomic level for primary tax determination
	 * @return Taxonomic ID (currently always returns -1 due to assertion)
	 */
	public int primaryTax(int level){
		//I have no idea how to implement this...
		IntHashMap map=new IntHashMap();
		assert(false);
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Print Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** String used to separate records in output formatting */
	private static String recordBreak="\n"; //"\n\n"
	
	/**
	 * Writes comparison results to a text stream using display parameters.
	 * Converts results to formatted text and outputs via the provided writer.
	 * @param params Display parameters controlling output format
	 * @param tsw Text stream writer for result output
	 */
	void writeResults(DisplayParams params, TextStreamWriter tsw){
		ByteBuilder sb=toText(params);
		tsw.print(sb);
	}
	
	/**
	 * Converts comparison results to formatted text representation.
	 * Handles SSU alignment if present, supports JSON output format,
	 * and applies various formatting modes including query-ref-ANI and
	 * constellation formats. Limits output to maximum record count.
	 *
	 * @param params Display parameters controlling output format and limits
	 * @return ByteBuilder containing formatted comparison results
	 */
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
	
	/**
	 * Initiates SSU (Small Subunit) sequence alignment for comparison results.
	 * Uses the alignment thread pool to process SSU alignments for up to
	 * maxRecords comparisons. Only processes if the query sketch contains SSU data.
	 * @param maxRecords Maximum number of comparisons to align
	 */
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
	
	/** The query sketch associated with these comparison results */
	public final Sketch sketch;
	/** List of reference sketches used in comparisons */
	public ArrayList<Sketch> refSketchList;
	/** Two-dimensional array tracking taxonomic hits across comparison levels */
	public int[][] taxHits;
	/** Sorted list of comparison results between query and reference sketches */
	public ArrayList<Comparison> list;
	/** Total number of comparison records processed */
	public int totalRecords=0;
	
}
