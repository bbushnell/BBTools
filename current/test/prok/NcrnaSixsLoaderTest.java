package test.prok;

import prok.*;

import java.io.File;
import java.util.Arrays;
import consensus.BaseGraph;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import parse.LineParser1;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Group binding, independent old-resource parity and resolver isolation checks.
 * @author Brian, Raiden */
public final class NcrnaSixsLoaderTest {
	public static void main(String[] args){
		String baseline=null, outdir=null;
		for(String arg:args){
			if(arg.startsWith("baseline=")){baseline=arg.substring(9);}
			else if(arg.startsWith("outdir=")){outdir=arg.substring(7);}
			else{throw new IllegalArgumentException("Expected baseline= and outdir=: "+arg);}
		}
		check(baseline!=null && outdir!=null && new File(outdir).isDirectory(), "Retained baseline and fresh evidence directory are required");
		Shared.setThreads(1);Read.TO_UPPER_CASE=true;
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final NcrnaFamily first=checkFamily(baseline, "sixs_rf00013", "RF00013", resolver);
		final String firstName=first.modelNames[0];
		checkFamily(baseline, "sixs_rf01685", "RF01685", resolver);
		check(first.modelNames[0].equals(firstName) && first.library.length==590, "Loading another member must preserve the first container");
		final String manifest=Data.findPath("?callgenes_6S/manifest.tsv", false);
		check(manifest!=null, "The paired package requires one parent manifest");
		final File root=new File(manifest).getParentFile();
		check(!new File(root, "RF00013/manifest.tsv").exists() && !new File(root, "RF01685/manifest.tsv").exists(), "Members must not acquire independent manifests");
		for(String mode:new String[]{"missing", "unknown", "duplicate", "escape", "boundary"}){
			final File altered=new File(outdir, mode);check(altered.mkdir(), "Use a fresh malformed-manifest fixture directory");
			mutate(manifest, new File(altered, "manifest.tsv").toString(), mode);
			final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();options.resources.root=altered.toString();
			final String reason=mode.equals("missing") ? "RF01685/models" : mode.equals("unknown") ? "Unknown group member"
				: mode.equals("duplicate") ? "Duplicate family resource role" : mode.equals("escape") ? "escapes" : "Unsupported resource role";
			fails(()->NcrnaFamilyLoader.load("sixs_rf00013", options, new FamilyResourceResolver()), reason);
		}
		final FamilyResourceResolver.Overrides wrongFamily=new FamilyResourceResolver.Overrides();
		final String r58Manifest=Data.findPath("?callgenes_5.8S/manifest.tsv", false);
		check(r58Manifest!=null, "The retained singleton package is required to test header isolation");
		wrongFamily.root=new File(r58Manifest).getParent();
		resolver.resolve("r58", "consensus", wrongFamily);
		fails(()->resolver.resolve("lsu", "consensus", wrongFamily), "Unexpected header:");
		System.out.println("PASS grouped6S: exact member resources/policies, one manifest, cross-member reuse, full-group schema failures and explicit-root identity isolation");
	}
	private static NcrnaFamily checkFamily(String baseline, String id, String member, FamilyResourceResolver resolver){
		final boolean rf13=id.equals("sixs_rf00013");
		final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.get(id);
		check(entry.group!=null && entry.member.equals(member), "Member identity must be explicit in registry bindings");
		final NcrnaFamily f=NcrnaFamilyLoader.load(id, new NcrnaFamilyLoader.Options(), resolver);
		check(resolver.resolve(id, "consensus", null).path.endsWith("/callgenes_6S/"+member+"/consensus.fa"), "Resolver must construct member paths from the group root");
		final String prefix=baseline+"/resources/"+id;
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(prefix+"_consensus.fa");
		final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		final BaseGraph[] models=TrnaConsensusBuilder.loadModels(prefix+"_models.hbmc.gz");
		check(Arrays.deepEquals(library, f.library) && Arrays.equals(names, f.modelNames) && models.length==f.models.length, "Ordered member libraries changed");
		for(int i=0; i<models.length; i++){check(models[i].toText().toString().equals(f.models[i].toText().toString()), "Decoded member HBM changed at "+i);}
		final String seeds=prefix+"_17mers.fa"+(rf13 ? ".gz" : "");
		keys(ProkObject.loadLongKmers(seeds, 17), f.kmerSet);keys(f.kmerSet, ProkObject.loadLongKmersByType(17, id));
		check(f.kLong==17 && f.seedMinHits==1 && !f.seedDistinct && f.minLen==60 && f.windowPad==(rf13 ? 250 : 125)
			&& f.indexK==7 && f.indexTopN==60 && f.fixedMinHits==1 && !f.adaptive && f.adaptFloor==0f && f.adaptTopFrac==0f && f.adaptQFrac==0f
			&& f.scoreA==0f && f.scoreB==20f && f.idPass==.70f && f.idBorderline==(rf13 ? .60f : .70f)
			&& f.hbmPass==(rf13 ? .64f : .75f) && f.collapseFrac==.85f && f.outputType==ProkObject.RNA
			&& f.scavengePass2 && f.trimAlignmentExtent && !f.reuseConsensusAlignment && f.boundary5NetTemplate==null,
			"Grouped layout must not alter either member's measured policy");
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.enabled=false;options.resources.root=baseline+"/absent-group";
		check(NcrnaFamilyLoader.load(id, options, resolver)==null, "Disabled members must not resolve the group");
		options.enabled=true;options.resources.root=null;options.boundary=true;
		options.resources.boundary=new File(baseline, "absent_boundary_"+System.nanoTime()+".bbnet").toString();
		fails(()->NcrnaFamilyLoader.load(id, options, resolver), "boundary");
		options.boundary=false;options.resources.seeds=seeds;options.windowPad=300;options.topN=2;options.minHits=3;
		final NcrnaFamily custom=NcrnaFamilyLoader.load(id, options, resolver);
		check(custom.windowPad==300 && custom.indexTopN==2 && custom.fixedMinHits==3, "Per-member overrides must remain independent");
		keys(f.kmerSet, custom.kmerSet);
		return f;
	}
	/** Malformed rows retain real provenance fields; only the deliberately tested binding changes. */
	private static void mutate(String source, String dest, String mode){
		check(!new File(dest).exists(), "Manifest fixtures must not overwrite previous evidence");
		final ByteFile in=ByteFile.makeByteFile(source, false);
		final ByteStreamWriter out=new ByteStreamWriter(dest, false, false, false);out.start();
		final LineParser1 p=new LineParser1('\t');
		try{for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);
			if(line.length==0 || line[0]=='#' || p.termEquals("member", 0)){out.print(new ByteBuilder().append(line).nl());continue;}
			if(mode.equals("missing") && p.termEquals("RF01685", 0) && p.termEquals("models", 1)){continue;}
			final boolean changed=p.termEquals("RF00013", 0) && p.termEquals("consensus", 1);
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<8; i++){
				if(i>0){row.append('\t');}
				if(changed && mode.equals("unknown") && i==0){row.append("UNKNOWN");}
				else if(changed && mode.equals("boundary") && i==1){row.append("unknown_boundary_role");}
				else{if(changed && mode.equals("escape") && i==4){row.append("../");}p.appendTerm(row, i);}
			}
			row.nl();out.print(row);
			if(changed && mode.equals("duplicate")){out.print(row);}
		}}finally{check(!in.close(), "Manifest fixture input read failed");check(!out.poisonAndWait(), "Manifest fixture write failed");}
	}
	private static void keys(LongHashSet a, LongHashSet b){
		final long[] x=a.toArray(), y=b.toArray();Arrays.sort(x);Arrays.sort(y);
		check(Arrays.equals(x, y), "Distinct subtype forward seeds must match the retained baseline");
	}
	private static void fails(Runnable action, String reason){
		try{action.run();}catch(IllegalArgumentException expected){check(expected.getMessage().contains(reason), "Wrong failure: "+expected);return;}
		throw new AssertionError("Expected group/resource failure: "+reason);
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
