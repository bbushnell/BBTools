package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import ddl.DDLFamilyProfiles;
import ddl.DDLPentamerBlacklist;
import ddl.ProteinKmerGeometry;

import prok.StopTailCascadeResources;

import prot.FamilyArtifactRebind;
import prot.FamilyCalibration.VerificationMode;
import prot.FamilyScoringPackage;
import prot.FamilyTailCascade.Assignment;

import shared.Shared;

import test.prok.StopTailCascadeStructuralTest;

/** Four-artifact inverse proof, actual tiny loaders and pre-reader rejection fixtures.
 * @author Sayu
 */
public final class FamilyArtifactRebindTest{
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh absolute fixture root");
		final Path root=Paths.get(args[0]);
		check(root.isAbsolute() && !Files.exists(root), "Fresh absolute fixture root required");
		Files.createDirectory(root);
		Shared.setThreads(1);
		check(java.nio.charset.Charset.defaultCharset().equals(StandardCharsets.UTF_8), "UTF-8 locale is required for the terminal filename fixture");
		final byte[] query=new byte[100];
		query[0]='M';
		final String acids="ACDEFGHIKLMNPQRSTVWY";
		for(int i=1; i<query.length; i++){
			query[i]=(byte)acids.charAt((i*7+i/20)%20);
		}
		final Path all=StopTailCascadeStructuralTest.writeInformationResources(root.resolve("resources"), query, false);
		final Path child=StopTailCascadeResources.readStructural(all.toString(), "unmeasured")
		.resolvedResourceBindings().get("short40.manifest");
		final Path longConfig=all.getParent().resolve("classifier.config");
		final Path shortConfig=child.getParent().resolve("classifier.config");
		final Path oldList=write(root.resolve("old.faa"), ">fixture\nAAAAA\n");
		final Path newList=write(root.resolve("new.faa"), ">fixture\nAAAAA\n");
		bindFixtureBlacklist(shortConfig, Arrays.copyOfRange(query, 1, query.length), oldList);
		bindFixtureBlacklist(longConfig, query, oldList);
		final Path[] sources={
			shortConfig, longConfig, child, all
		};
		for(int i=0; i<sources.length; i++){
			final String delimiter=i<2?"=":"\t";
			String text=read(sources[i])+"# Preserve UTF-8 comment: \u00e9\nblacklist"+delimiter+oldList+"\nblacklistsha80"+delimiter+"unmeasured\n"+
				"blacklistk"+delimiter+"5\nblacklistmode"+delimiter+"reduced_literal\nblacklistwords"+delimiter+"1\nblacklistprovenance"+delimiter+"fixture\n";
			if(i==1){
				text=text.substring(0, text.length()-1).replace("\n", "\r\n");
			}
			write(sources[i], text);
		}
		final String[] targets={
				"short.classifier.config", "long.classifier.config", "short40.tsv", "all.tsv"
		};
		final byte[][] originals=new byte[4][];
		for(int i=0; i<4; i++){
			originals[i]=Files.readAllBytes(sources[i]);
		}
		final Path out=root.resolve("rebound");
		final String batch=batch(sources, targets, oldList, newList, child, out);
		final Path table=write(root.resolve("bindings.tsv"), batch);
		FamilyArtifactRebind.rebind(table, out);
		check(read(out.resolve("PASS")).equals("FAMILY_ARTIFACT_REBIND_PASS artifacts=4 changed_path_fields=5 classifier_runs=0\n"), "Exact four/five completion");
		for(int i=0; i<4; i++){
			final byte[] current=Files.readAllBytes(out.resolve(targets[i]));
			String inverse=new String(current, StandardCharsets.UTF_8).replace(newList.toString(), oldList.toString());
			if(i==3){
				inverse=inverse.replace(out.resolve("short40.tsv").toString(), child.toString());
			}
			check(Arrays.equals(originals[i], inverse.getBytes(StandardCharsets.UTF_8)), "Whole-file inverse with unchanged CRLF and terminal state: "+i);
			check(Arrays.equals(originals[i], Files.readAllBytes(sources[i])), "Original artifact untouched: "+i);
		}
		final StopTailCascadeResources before=StopTailCascadeResources.readStructural(all.toString(), "unmeasured");
		final StopTailCascadeResources after=StopTailCascadeResources.readStructural(out.resolve("all.tsv").toString(), "unmeasured");
		check(before.firstStart==after.firstStart && before.allFamilySets==after.allFamilySets && before.marginScale==after.marginScale,
			"Caller policy remains unchanged");
		final Map<String, Path> detached=after.resolvedResourceBindings();
		detached.clear();
		check(!after.resolvedResourceBindings().isEmpty(), "Role view is detached");
		same(before.load(1).newWorker().assign("before", query), after.load(1).newWorker().assign("after", query));
		for(int i=0; i<2; i++){
			final FamilyScoringPackage a=new FamilyScoringPackage(sources[i], VerificationMode.STRUCTURAL);
			final FamilyScoringPackage b=new FamilyScoringPackage(out.resolve(targets[i]), VerificationMode.STRUCTURAL);
			final byte[] packageQuery=i==0?Arrays.copyOfRange(query, 1, query.length):query;
			same(a.load(1).newWorker().assign("before", packageQuery), b.load(1).newWorker().assign("after", packageQuery));
		}
		final Path unicodeOld=root.resolve("old-unicode-\u00e9");
		final Path unicodeNew=root.resolve("new-unicode-\u00e9");
		Files.copy(oldList, unicodeOld);
		Files.copy(oldList, unicodeNew);
		String unicodeText=read(longConfig).replace(oldList.toString(), unicodeOld.toString());
		unicodeText=unicodeText.replace("blacklist="+unicodeOld+"\r\n", "blacklist="+unicodeOld+" \t\r\n");
		final Path unicodeSource=write(root.resolve("unicode.config"), unicodeText);
		final Path unicodeOut=root.resolve("unicode-out");
		final Path unicodeTable=write(root.resolve("unicode.tsv"), header()+row("u", "classifier", unicodeSource, "u.config", "blacklist", unicodeOld.toString(), unicodeNew.toString()));
		FamilyArtifactRebind.rebind(unicodeTable, unicodeOut);
		check(new FamilyScoringPackage(unicodeOut.resolve("u.config"), VerificationMode.STRUCTURAL).resolvedOptions().get("blacklist").equals(unicodeNew.toString()),
			"Terminal UTF-8 filename retains its actual binding and classifier whitespace is trimmed");
		check(Arrays.equals(Files.readAllBytes(unicodeSource), read(unicodeOut.resolve("u.config")).replace(unicodeNew.toString(), unicodeOld.toString()).getBytes(StandardCharsets.UTF_8)),
			"Terminal UTF-8 path, CRLF and surrounding classifier whitespace survive whole-file inverse");
		check(read(unicodeOut.resolve("bindings.tsv")).contains(unicodeNew.toString()), "Receipt retains UTF-8 replacement bytes");
		final Path oldSpace=root.resolve("old.faa ");
		final Path newSpace=root.resolve("new.faa ");
		Files.copy(oldList, oldSpace);
		Files.copy(oldList, newSpace);
		final Path spaceSource=write(all.getParent().resolve("space-all.tsv"), read(all).replace("blacklist\t"+oldList+"\n", "blacklist\t"+oldSpace+"\n"));
		final Path spaceOut=root.resolve("space-out");
		final Path spaceTable=write(root.resolve("space.tsv"), header()+row("space", "cascade", spaceSource, "all.tsv", "blacklist", oldSpace.toString(), newSpace.toString()));
		FamilyArtifactRebind.rebind(spaceTable, spaceOut);
		check(Files.isSameFile(StopTailCascadeResources.readStructural(spaceOut.resolve("all.tsv").toString(), "unmeasured").resolvedResourceBindings().get("blacklist"), newSpace),
			"Cascade value spaces remain filename data");
		check(Arrays.equals(Files.readAllBytes(spaceSource), read(spaceOut.resolve("all.tsv")).replace(newSpace.toString(), oldSpace.toString()).getBytes(StandardCharsets.UTF_8)),
			"Cascade trailing-space path has exact inverse bytes");
		// Metadata-only GC fixtures deliberately contain no scientific rows. Graph
		// binding must discover/protect them without attempting a statistics load.
		final Path stats=write(all.getParent().resolve("gc.tsv"), "not loaded by metadata binding\n");
		final Path gcMap=write(all.getParent().resolve("gc-map.tsv"), "not loaded by metadata binding\n");
		final Path gcAll=write(all.getParent().resolve("gc-all.tsv"), read(all)+"familygc\t"+stats+"\nfamilygcsha80\tunmeasured\n"+
			"familygcmap\t"+gcMap+"\nfamilygcmapsha80\tunmeasured\nfamilygcfamilies\t61\nfamilygcmembers\t1830\n"+
			"familygclambda\t24\nfamilygcfilter\tt\nfamilygcz\t8\n");
		final Path gcOut=root.resolve("gc-out");
		final Path gcTable=write(root.resolve("gc-bindings.tsv"), header()+
		row("gc", "cascade", gcAll, "all.tsv", "blacklist", oldList.toString(), newList.toString()));
		FamilyArtifactRebind.rebind(gcTable, gcOut);
		final Map<String, Path> gcRoles=StopTailCascadeResources.readStructural(gcOut.resolve("all.tsv").toString(), "unmeasured").resolvedResourceBindings();
		check(Files.isSameFile(stats, gcRoles.get("familygc")) && Files.isSameFile(gcMap, gcRoles.get("familygcmap")), "GC dependencies retain their resolved roles");
		check(gcRoles.containsKey("short40.familyreferences") && gcRoles.containsKey("short40.dependency.0"), "Child profile dependencies are included");
		write(gcAll, read(gcAll).replace("familygc\t"+stats, "familygc\tgc.tsv"));
		final Path badGcOut=root.resolve("relative-gc-out");
		final Path badGc=write(root.resolve("relative-gc.tsv"), header()+
		row("gc", "cascade", gcAll, "all.tsv", "blacklist", oldList.toString(), newList.toString())+
		row("other", "classifier", shortConfig, "gc.tsv", "blacklist", oldList.toString(), newList.toString()));
		rejectMessage(()->FamilyArtifactRebind.rebind(badGc, badGcOut), badGcOut, "familygc");
		reject(()->FamilyArtifactRebind.rebind(table, out), out, "Existing output root");
		bad(root, "wrong-old", batch.replace(oldList.toString(), root.resolve("wrong.faa").toString()));
		bad(root, "duplicate-target", batch.replace("long.classifier.config", "short.classifier.config"));
		bad(root, "unsafe-target", batch.replace("long.classifier.config", "../escape"));
		bad(root, "scientific-key", batch.replace("\tblacklist\t", "\ttop\t"));
		bad(root, "duplicate-key", batch+batch.substring(batch.indexOf('\n', batch.indexOf('\n')+1)+1).split("\n")[0]+"\n");
		// A complete, valid-looking destination exists, but relocation must not silently
		// reinterpret the LONG reference as the SHORT classifier's output file.
		final String prior=read(longConfig);
		write(longConfig, prior.replace(all.getParent().resolve("references.tsv").toString(), "references.tsv"));
		final Path relativeOut=root.resolve("relative-out");
		final String relative=batch(sources, new String[]{
				"references.tsv", targets[1], targets[2], targets[3]
		}, oldList, newList, child, relativeOut);
		reject(()->FamilyArtifactRebind.rebind(write(root.resolve("relative.tsv"), relative), relativeOut), relativeOut, "Resolved relative binding changed");
		write(longConfig, prior);
		// These inputs also lack mandatory classifier fields: the specific rejection
		// proves no ordinary reader was entered (or diagnostic matrix hashed) first.
		for(String field:new String[]{
				"uncalibratedmatrixfile", "familymeans"
		}){
			final Path unsupported=write(root.resolve(field+".config"), "blacklist="+oldList+"\n"+field+"=/missing/payload\n");
			final Path rejectedOut=root.resolve(field+"-out");
			final Path input=write(root.resolve(field+".tsv"), header()+row("x", "classifier", unsupported, "x.config", "blacklist", oldList.toString(), newList.toString()));
			rejectMessage(()->FamilyArtifactRebind.rebind(input, rejectedOut), rejectedOut, "before structural reader");
			check(!Files.exists(rejectedOut), "Pre-reader failure creates no output root");
		}
		final Path fitted=write(root.resolve("fitted.tsv"), "#bbgap_tail_fitted\t1\nblacklist\t"+oldList+"\n");
		final Path fittedOut=root.resolve("fitted-out");
		final Path fittedTable=write(root.resolve("fitted-bindings.tsv"), header()+
		row("fitted", "cascade", fitted, "x.tsv", "blacklist", oldList.toString(), newList.toString()));
		rejectMessage(()->FamilyArtifactRebind.rebind(fittedTable, fittedOut), fittedOut, "Unsupported cascade header");
		final Path cycle=root.resolve("cycle.manifest.tsv");
		write(cycle, "#bbgap_tail_cascade\t7\nshort40manifest\t"+cycle+"\nblacklist\t"+oldList+"\n");
		final Path cycleOut=root.resolve("cycle-out");
		final Path cycleTable=write(root.resolve("cycle.tsv"), header()+row("x", "cascade", cycle, "x.tsv", "blacklist", oldList.toString(), newList.toString()));
		rejectMessage(()->FamilyArtifactRebind.rebind(cycleTable, cycleOut), cycleOut, "Cyclic");
		final Path orphanOut=root.resolve("orphan-out");
		final Path orphan=write(root.resolve("orphan.tsv"), header()+row("s", "cascade", child, "s.tsv", "blacklist", oldList.toString(), newList.toString()));
		reject(()->FamilyArtifactRebind.rebind(orphan, orphanOut), orphanOut, "SHORT needs composite context");
		System.out.println("FAMILY_ARTIFACT_REBIND_TEST_PASS checks="+checks+" artifacts=4 edits=5");
	}

	/** Bind the inert fixture blacklist only after proving every saved modal row still matches it. */
	private static void bindFixtureBlacklist(Path config, byte[] member, Path blacklistPath)throws Exception{
		final Map<String, String> options=new FamilyScoringPackage(config, VerificationMode.STRUCTURAL).resolvedOptions();
		final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(Integer.parseInt(options.get("k")), options.get("aabits"));
		final DDLPentamerBlacklist blacklist=DDLPentamerBlacklist.loadReducedStructural(blacklistPath, "unmeasured", 1, geometry, "fixture");
		final char[] unfiltered=new char[512], filtered=new char[512];
		DDLFamilyProfiles.sketchPositional(member, geometry, 0, 4, 2, unfiltered, null);
		DDLFamilyProfiles.sketchPositional(member, geometry, 0, 4, 2, filtered, blacklist);
		check(Arrays.equals(unfiltered, filtered), "Fixture blacklist leaves the actual full member sketch unchanged");
		int retained=0;
		final StringBuilder modal=new StringBuilder();
		for(int i=0; i<filtered.length; i++){
			if(filtered[i]!=0 && ++retained>3){
				filtered[i]=0;
			}
			if(i>0){
				modal.append(',');
			}
			modal.append((int)filtered[i]);
		}
		check(retained>=3, "The existing fixture retains three registers");
		final List<String> list=Files.readAllLines(Paths.get(options.get("ddlprofiles")), StandardCharsets.UTF_8);
		check(list.size()==1, "One literal fixture profile shard");
		final String[] entry=list.get(0).split("\t", -1);
		check(entry.length==3, "Profile and summary identities remain explicit");
		final Path profile=Paths.get(entry[0]);
		final List<String> rows=Files.readAllLines(profile, StandardCharsets.UTF_8);
		check(rows.size()==62 && rows.get(0).equals(DDLFamilyProfiles.HEADER16), "The original61-family fixture profile schema remains intact");
		for(int i=1; i<rows.size(); i++){
			final String[] fields=rows.get(i).split("\t", -1);
			check(fields.length==10 && fields[7].equals(modal.toString()), "Every retained modal row equals the blacklist-aware fixture sketch");
		}
		final Path summary=profile.getParent().resolve("summary.tsv");
		final List<String> lines=Files.readAllLines(summary, StandardCharsets.UTF_8);
		check(lines.size()==2 && !lines.get(0).contains("blacklist_"), "Bind only the original unfiltered fixture summary");
		write(summary, lines.get(0)+"\tblacklist_sha80\tblacklist_words\t"+blacklist.metadataColumns()+"\n"+
			lines.get(1)+'\t'+blacklist.sha80()+'\t'+blacklist.size()+'\t'+blacklist.metadataValues()+"\n");
	}

	/** Describe exactly four independent blacklist edits and the parent-to-child relocation. */
	private static String batch(Path[] sources, String[] targets, Path oldList, Path newList, Path child, Path out){
		String text=header();
		for(int i=0; i<4; i++){
			text+=row("a"+i, i<2?"classifier":"cascade", sources[i], targets[i], "blacklist", oldList.toString(), newList.toString());
		}
		return text+row("a3", "cascade", sources[3], targets[3], "short40manifest", child.toString(), out.resolve("short40.tsv").toString());
	}

	private static String header(){
		return "#bbgap_family_rebind\t1\nartifact\tformat\tsource\ttarget\tkey\texpected\treplacement\n";
	}

	private static String row(String id, String format, Path source, String target, String key, String oldValue, String replacement){
		return id+'\t'+format+'\t'+source+'\t'+target+'\t'+key+'\t'+oldValue+'\t'+replacement+'\n';
	}

	/** Compare actual tiny-loader decisions and work counters, not just parsed path maps. */
	private static void same(Assignment a, Assignment b){
		check(a.assigned()==b.assigned() && a.terminal.equals(b.terminal) && a.rawHitFloor==b.rawHitFloor && a.ddlMatches==b.ddlMatches &&
		a.dpAttempts==b.dpAttempts && a.admittedFamilies==b.admittedFamilies && Float.floatToIntBits(a.identity)==Float.floatToIntBits(b.identity) &&
		Float.floatToIntBits(a.blosum)==Float.floatToIntBits(b.blosum), "Actual tiny classifier assignments/counters match");
	}

	private interface Action{
		void run()throws Exception;
	}

	private static void bad(Path root, String name, String text)throws Exception{
		final Path table=write(root.resolve(name+".tsv"), text);
		final Path out=root.resolve(name+"-out");
		reject(()->FamilyArtifactRebind.rebind(table, out), out, name);
	}

	private static void reject(Action action, Path out, String why)throws Exception{
		rejectMessage(action, out, null);
	}

	private static void rejectMessage(Action action, Path out, String fragment)throws Exception{
		final boolean exists=Files.exists(out.resolve("PASS"));
		try{
			action.run();
		}catch(IllegalArgumentException | IllegalStateException e){
			check(fragment==null || e.getMessage()!=null && e.getMessage().contains(fragment), "Wrong rejection: "+e);
			check(exists || !Files.exists(out.resolve("PASS")), "Failed operation published PASS");
			return;
		}
		throw new AssertionError("Expected rejection");
	}

	private static Path write(Path file, String text)throws Exception{
		return Files.write(file, text.getBytes(StandardCharsets.UTF_8));
	}

	private static String read(Path file)throws Exception{
		return new String(Files.readAllBytes(file), StandardCharsets.UTF_8);
	}

	private static void check(boolean valid, String why){
		if(!valid){
			throw new AssertionError(why);
		}
		checks++;
	}

	private static int checks;
}
