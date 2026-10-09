package ml;

import java.util.Arrays;

import json.JsonObject;
import json.JsonParser;
import structures.ByteBuilder;

/** Preserves output ownership when ModifyNN appends hidden nodes and heads.
 * @author Nilou
 */
final class ModifyNNPartition {

	static ModifyNNPartition grow(CellNet parent, int[] dims, int privatePerHead){
		final String text=parent.getTag("output_partition");
		if(text==null){
			if(privatePerHead!=0){throw new IllegalArgumentException("privateperhead requires output_partition metadata");}
			return null;
		}
		if(parent.dims.length<3){throw new IllegalArgumentException("Output partition requires a hidden layer");}
		final String trimmed=text.trim();
		if(!new JsonParser(trimmed).validate()){throw new IllegalArgumentException("Malformed output_partition JSON");}
		final JsonObject json=new JsonParser(trimmed).parseJsonObject();
		if(json==null || json.omap==null){throw new IllegalArgumentException("Output partition must be a JSON object");}
		final int oldWidth=parent.dims[parent.dims.length-2], oldHeads=parent.dims[parent.dims.length-1];
		if(integer(json.omap.get("hidden_width"))!=oldWidth || integer(json.omap.get("outputs"))!=oldHeads){
			throw new IllegalArgumentException("Output partition dimensions disagree with parent network");
		}
		final boolean[] oldShared=new boolean[oldWidth];
		final boolean[][] oldAllowed=new boolean[oldHeads][oldWidth];
		final Object scheme=json.omap.get("scheme");
		if("half_shared_v1".equals(scheme)){
			final int shared=integer(json.omap.get("shared")), privateCount=oldWidth-shared;
			if(shared!=oldWidth/2 || privateCount<oldHeads){throw new IllegalArgumentException("Invalid half-shared dimensions");}
			final Object[] counts=array(json.omap.get("private_counts")), order=array(json.omap.get("head_order"));
			if(counts.length!=oldHeads || order.length!=oldHeads){throw new IllegalArgumentException("Missing private counts or head order");}
			Arrays.fill(oldShared, 0, shared, true);
			int start=shared;
			for(int h=0; h<oldHeads; h++){
				final int count=integer(counts[h]), expected=privateCount/oldHeads+(h<privateCount%oldHeads ? 1 : 0);
				if(count!=expected || integer(order[h])!=h){throw new IllegalArgumentException("Half-shared head layout is inconsistent: head="+h);}
				Arrays.fill(oldAllowed[h], 0, shared, true);
				Arrays.fill(oldAllowed[h], start, start+count, true);
				start+=count;
			}
		}else if("modified_explicit_v1".equals(scheme)){
			readIndexes(json.omap.get("shared_inputs"), oldShared);
			final Object[] rows=array(json.omap.get("allowed_inputs"));
			if(rows.length!=oldHeads){throw new IllegalArgumentException("Output partition must declare every head");}
			for(int h=0; h<oldHeads; h++){
				readIndexes(rows[h], oldAllowed[h]);
				for(int i=0; i<oldWidth; i++){
					if(oldShared[i] && !oldAllowed[h][i]){throw new IllegalArgumentException("Declared shared node absent from head "+h);}
				}
			}
		}else{throw new IllegalArgumentException("Unknown output partition scheme: "+scheme);}
		validateEdges(parent, oldAllowed);
		final int width=dims[dims.length-2], heads=dims[dims.length-1], addedHeads=heads-oldHeads;
		if(privatePerHead>0 && (addedHeads<1 || (long)width-oldWidth!=(long)addedHeads*privatePerHead)){
			throw new IllegalArgumentException("Private growth requires exactly privateperhead appended nodes per new head");
		}
		final boolean[] shared=Arrays.copyOf(oldShared, width);
		final boolean[][] allowed=new boolean[heads][width];
		for(int h=0; h<oldHeads; h++){System.arraycopy(oldAllowed[h], 0, allowed[h], 0, oldWidth);}
		if(privatePerHead==0){
			Arrays.fill(shared, oldWidth, width, true);
			for(int h=0; h<oldHeads; h++){Arrays.fill(allowed[h], oldWidth, width, true);}
		}
		for(int h=oldHeads; h<heads; h++){
			System.arraycopy(shared, 0, allowed[h], 0, width);
			if(privatePerHead>0){
				final int start=oldWidth+(h-oldHeads)*privatePerHead;
				Arrays.fill(allowed[h], start, start+privatePerHead, true);
			}
		}
		return new ModifyNNPartition(shared, allowed);
	}

	private ModifyNNPartition(boolean[] shared_, boolean[][] allowed_){
		assert(allowed_.length>0 && shared_.length>0) : "Validated output partition must cover heads and last-hidden nodes";
		shared=shared_; allowed=allowed_;
	}

	/** Rejects forbidden active edges before growth can conceal invalid parent ownership. */
	private static void validateEdges(CellNet parent, boolean[][] allowed){
		final Cell[] layer=parent.net[parent.net.length-1];
		assert(layer.length==allowed.length) : "Partition head count was checked against parent dimensions";
		for(int h=0; h<layer.length; h++){
			final Cell cell=layer[h];
			if(cell.inputs==null){
				for(int i=0; i<cell.weights.length; i++){
					if(!allowed[h][i] && (cell.weights[i]!=0 || !parent.zeroAbsent())){throw new IllegalArgumentException("Forbidden partition weight at head="+h+" input="+i);}
				}
			}else{
				for(int i : cell.inputs){
					if(!allowed[h][i]){throw new IllegalArgumentException("Forbidden stored partition edge at head="+h+" input="+i);}
				}
			}
		}
	}

	private static Object[] array(Object value){
		if(!(value instanceof Object[])){throw new IllegalArgumentException("Partition index/count list must be an array");}
		return (Object[])value;
	}

	private static int integer(Object value){
		if(!(value instanceof Long)){throw new IllegalArgumentException("Partition dimensions and indexes must be JSON integers: "+value);}
		final long n=((Long)value).longValue();
		if(n<0 || n>Integer.MAX_VALUE){throw new IllegalArgumentException("Partition integer outside supported range: "+n);}
		return (int)n;
	}

	private static void readIndexes(Object value, boolean[] membership){
		int previous=-1;
		for(Object item : array(value)){
			final int i=integer(item);
			if(i<=previous || i>=membership.length){throw new IllegalArgumentException("Partition indexes must be sorted, unique, and in range");}
			membership[i]=true; previous=i;
		}
	}

	/** One-line metadata matches GPU modified_explicit_v1 ownership semantics. */
	String metadata(){
		final ByteBuilder out=new ByteBuilder();
		out.append("{\"allowed_inputs\":[");
		for(int h=0; h<allowed.length; h++){
			if(h>0){out.append(',');}
			appendIndexes(out, allowed[h]);
		}
		out.append("],\"hidden_width\":").append(shared.length).append(",\"outputs\":").append(allowed.length);
		out.append(",\"scheme\":\"modified_explicit_v1\",\"shared_inputs\":");
		appendIndexes(out, shared);
		return out.append('}').toString();
	}

	private static void appendIndexes(ByteBuilder out, boolean[] membership){
		assert(out!=null) : "Partition metadata appends indexes into one shared output buffer";
		out.append('['); boolean comma=false;
		for(int i=0; i<membership.length; i++){
			if(membership[i]){
				if(comma){out.append(',');}
				out.append(i); comma=true;
			}
		}
		out.append(']');
	}

	final boolean[] shared;
	final boolean[][] allowed;
}
