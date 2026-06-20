package bin;

import java.io.File;
import java.util.ArrayList;

import fileIO.ReadWrite;
import json.JsonObject;
import parse.PreParser;
import sketch.DisplayParams;
import sketch.SendSketch;
import sketch.Sketch;
import sketch.SketchIdComparator;
import sketch.SketchObject;
import sketch.SketchTool;

/**
 * Renames files with top sketch hit taxid.
 * Loads sequence files as sketches, compares them against RefSeq database,
 * and renames each file by prepending the taxonomic ID of the best match.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class FileRenamer {

	public static void main(String[] args) {
		PreParser pp=new PreParser(args, null, false);
		args=pp.args;
		DisplayParams params=new DisplayParams();
		params.format=DisplayParams.FORMAT_JSON;
		SketchObject.postParse();
		SketchTool tool=new SketchTool(10000, params);
		ArrayList<Sketch> inSketches=tool.loadSketches_MT(params, args);
		assert(inSketches.size()==args.length) : inSketches.size()+" != "+args.length;
		final int numLoaded=(inSketches.size());
		if(numLoaded>1 && params.mode==Sketch.PER_FILE){
			inSketches.sort(SketchIdComparator.comparator);//Otherwise they come out disordered
		}
		
		ArrayList<JsonObject> results=SendSketch.sendSketches(inSketches, "refseq", params);
		assert(results.size()==inSketches.size()) : inSketches.size()+" != "+results.size();
		for(int i=0; i<args.length; i++) {
			String fname=args[i];
			JsonObject result=results.get(i);
			
			JsonObject top=null;
			if(result!=null && result.jmapSize()>0) {
				for(String key : result.jmap.keySet()){
					top=result.jmap.get(key);
					break;
				}
			}
			SketchRecord topHit=(top==null ? null : new SketchRecord(top));
			int taxid=(topHit==null ? -1 : topHit.taxid);
			
			File f=new File(fname);
			//[bin/FileRenamer#001 FIXED 2026-06-20]: prepend tid to the BASENAME, not the whole path. Was "tid_X_"+fname,
			//which for a path arg (e.g. data/x.fa) produced "tid_X_data/x.fa" (a bogus leading dir) -> renameTo silently failed.
			String dir=ReadWrite.getPath(fname);
			String fname2=(dir==null ? "" : dir)+"tid_"+taxid+"_"+ReadWrite.stripPath(fname);
			File f2=new File(fname2);
			assert(f.exists());
			assert(!f2.exists());
			boolean renamed=f.renameTo(f2);//[bin/FileRenamer#001] was: return value ignored -> silent failure; crash-loud now
			assert(renamed) : "Failed to rename "+fname+" -> "+fname2;
		}
	}
	
}
