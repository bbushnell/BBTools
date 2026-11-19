package test.benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import benchmark.CallGenesReplaySelectionCache;
import benchmark.CallGenesReplaySelectionCache.Inventory;
import structures.ByteBuilder;

/** Relocation metadata regressions with synthetic payloads and retained strict controls.
 * Does not call genomes or train models. Synthetic identities are inherited labels;
 * the ordinary strict-reader controls measure existing manifest hashes, but do not
 * create or rewrite their identity labels. The execution wrapper must separately check
 * payload byte identity; a deliberately wrong same-size payload is written for that check.
 * Run serially in a fresh output directory on a filesystem supporting symlinks.
 * Generated files are retained for inspection, including mutated inputs if a test aborts.
 * @author Qiqi, Fischl
 */
public final class CallGenesReplayModelRelocationTest{
	/*--------------------------------------------------------------*/
	/*----------------         Entry point          ----------------*/
	/*--------------------------------------------------------------*/

	/** Require -ea and explicit out/inventory/cache/inventorypin/cachepin arguments.
	 * Create an absent absolute output directory, run synthetic cases, then check the
	 * retained 24-model/983-genome strict manifests and altered-pin rejections.
	 * Write and print PASS only after every check; input manifests remain read-only. */
	public static void main(String[] args)throws Exception{
		boolean assertions=false;
		assert(assertions=true);
		check(assertions, "Run relocation fixtures with -ea");
		Path root=null, originalInventory=null, originalCache=null;
		String inventoryPin=null, cachePin=null;
		final java.util.HashSet<String> seen=new java.util.HashSet<String>();
		for(String arg:args){
			final int split=arg.indexOf('=');
			check(split>0 && split<arg.length()-1, "Require nonempty key=value arguments");
			final String key=arg.substring(0, split), value=arg.substring(split+1);
			check(seen.add(key), "Duplicate fixture argument");
			if(key.equals("out")){
				root=Paths.get(value);
			}else if(key.equals("inventory")){
				originalInventory=Paths.get(value);
			}else if(key.equals("cache")){
				originalCache=Paths.get(value);
			}else if(key.equals("inventorypin")){
				inventoryPin=value;
			}else if(key.equals("cachepin")){
				cachePin=value;
			}else{
				throw new IllegalArgumentException("Unknown fixture argument: "+key);
			}
		}
		check(root!=null && root.isAbsolute() && originalInventory!=null && originalCache!=null
			&& inventoryPin!=null && cachePin!=null,
			"Explicit fixture output and retained strict-control manifests are required");
		Files.createDirectory(root);
		new CallGenesReplayModelRelocationTest(root).run();
		// These unchanged APIs hash retained manifests against their existing labels.
		// The fixture neither generates replacement labels nor rewrites those files.
		final Inventory strict=CallGenesReplaySelectionCache.readInventory(originalInventory, inventoryPin);
		check(!strict.relocated(), "Original reader unexpectedly enabled relocation");
		final CallGenesReplaySelectionCache cache=CallGenesReplaySelectionCache.read(originalCache, cachePin, strict);
		check(strict.size()==24 && cache.size()==983, "Retained strict control must cover24 models and983 genomes");
		final Path inv=originalInventory, source=originalCache;
		final String badInventory=wrongPin(inventoryPin), badCache=wrongPin(cachePin);
		reject(()->CallGenesReplaySelectionCache.readInventory(inv, badInventory), "Manifest sha80 mismatch");
		reject(()->CallGenesReplaySelectionCache.read(source, badCache, strict), "Manifest sha80 mismatch");
		Files.write(root.resolve("PASS"), ("MODEL_RELOCATION_METADATA_PASS checks="+checks+
			" strict_models=24 strict_genomes=983 caller_runs=0 payload_byte_identity=EXTERNAL_WORKER\n").getBytes(StandardCharsets.US_ASCII));
		System.out.print(new String(Files.readAllBytes(root.resolve("PASS")), StandardCharsets.US_ASCII));
	}

	/*--------------------------------------------------------------*/
	/*----------------       Synthetic setup        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create owned original/resource files and two-row inventory, relocation and selection metadata.
	 * The FASTA paths are names only; the three-byte payloads are not valid trained models. */
	private CallGenesReplayModelRelocationTest(Path root)throws Exception{
		assert(root!=null) : "Every generated fixture file belongs under the fresh output directory";
		resources=Files.createDirectory(root.resolve("resources"));
		final Path originals=Files.createDirectory(root.resolve("originals"));
		oldA=originals.resolve("a.pgm");
		oldB=originals.resolve("b.pgm");
		newA=resources.resolve("a.pgm");
		newB=resources.resolve("b.pgm");
		write(oldA, "ABC");
		write(oldB, "DEF");
		write(newA, "ABC");
		write(newB, "DEF");
		write(root.resolve("wrong_same_size.pgm"), "XYZ");
		inventory=root.resolve("inventory.tsv");
		mapping=root.resolve("map.tsv");
		cacheFile=root.resolve("cache.tsv");
		fasta=root.resolve("fixture.fa").toString();
		inventoryText="path\tsha80\tsize_bytes\n"+oldA+"\t"+LABEL+"\t3\n"+oldB+"\t"+LABEL+"\t3\n";
		mapText=MAP_HEADER+oldA+"\ta.pgm\n"+oldB+"\tb.pgm\n";
		cacheText="#CallGenesReplaySelectionCache\t1\tmodel_inventory_sha80\t"+LABEL+"\n"
			+"genome\tfasta\tfasta_sha80\tpolicy\tdomain\tphylum\tsource\tfallback\tcache_hit\tmodel_path\tmodel_sha80\n"
			+"g0000\t"+fasta+"\t"+LABEL+"\tPER_PHYLUM\tBacteria\tFixture\tserver\tnone\t0\t"+oldA+"\t"+LABEL+"\n"
			+"g0001\t"+fasta+".missing\t"+LABEL+"\tPER_PHYLUM\tBacteria\tMissing\tserver\tphylum_model_missing\t0\t"+oldB+"\t"+LABEL+"\n";
		write(inventory, inventoryText);
		write(mapping, mapText);
		write(cacheFile, cacheText);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Relocation checks       ----------------*/
	/*--------------------------------------------------------------*/

	/** Exercise explicit relocation, inverse bytes, roster/taxonomy matching and malformed maps.
	 * Expected negative cases restore synthetic metadata before continuing; unexpected
	 * failures abort and may leave that synthetic mutation in place. Accepted retained
	 * resources are never edited. No FASTA or model parsing is performed here. */
	private void run()throws Exception{
		Inventory models=load();
		check(models.size()==2 && models.resolvedPath(oldA.toString()).equals(newA.toString()), "Valid map resolution");
		models.requireFallback(newB.toString());
		reject(()->load().requireFallback(oldB.toString()), "Fallback is absent");
		CallGenesReplaySelectionCache cache=CallGenesReplaySelectionCache.readRelocatedCache(cacheFile, LABEL, models);
		check(cache.size()==2, "Complete fixture cache");
		cache.requireGenome(0, "g0000", fasta, LABEL);
		final CallGenesReplaySelectionCache valid=cache;
		reject(()->valid.requireGenome(0, "g0001", fasta, LABEL), "roster/FASTA identity changed");
		reject(()->valid.validateResolved(0, oldA.toString(), "Bacteria", "Fixture", "server", "none", false), "metadata differ");
		reject(()->valid.validateResolved(1, newB.toString(), "Bacteria", "Missing", "server", "none", false), "metadata differ");
		reject(()->valid.validateResolved(0, newA.toString(), "Bacteria", "Fixture", "server", "none", true), "metadata differ");
		ByteBuilder inverse=CallGenesReplaySelectionCache.begin(models),
			resolution=CallGenesReplaySelectionCache.beginResolution();
		cache.appendResolved(inverse, resolution, 0, newA.toString(), "Bacteria", "Fixture", "server", "none", false);
		cache.appendResolved(inverse, resolution, 1, newB.toString(), "Bacteria", "Missing", "server",
			"phylum_model_missing", false);
		cache.verifyInverse(inverse);
		check(new String(resolution.toBytes(), StandardCharsets.UTF_8).contains(newA.toString()),
			"Resolution evidence omits actual path");
		reject(()->valid.verifyInverse(new ByteBuilder("wrong\n")), "inverse differs");
		badMap(mapText.replace("Relocation\t1", "Relocation\t2"), "version header");
		badMap(MAP_HEADER+oldA+"\ta.pgm\n", "missing inventory");
		badMap(mapText+oldA+"\ta.pgm\n", "extra row");
		badMap(mapText.replace(oldB.toString(), oldA.toString()), "inventory order");
		badMap(mapText.replace("\tb.pgm\n", "\ta.pgm\n"), "Duplicate relocation");
		badMap(mapText.replace("\tb.pgm\n", "\t"+newB+"\n"), "normalized relative path");
		badMap(mapText.replace("\tb.pgm\n", "\t../outside.pgm\n"), "escapes resource root");
		badMap(mapText.replace("\tb.pgm\n", "\tmissing.pgm\n"), "Missing or aliased");
		badMap(mapText+"#CallGenesReplayModelRelocation\t1\n", "extra row");
		badMap(mapText.substring(0, mapText.length()-1), "terminal LF");
		final Path link=resources.resolve("alias.pgm");
		Files.createSymbolicLink(link, newB);
		badMap(mapText.replace("\tb.pgm\n", "\talias.pgm\n"), "contains a symlink");
		write(newB, "LONG");
		reject(()->load(), "size differs");
		write(newB, "DEF");
		badCache(cacheText.replace("model_inventory_sha80\t"+LABEL, "model_inventory_sha80\t"+OTHER),
			"model-set binding changed");
		badCache(cacheText.replace("\t"+oldA+"\t"+LABEL, "\t"+oldA+"\t"+OTHER), "model identity differs");
		badCache(cacheText.replace("\t0\t", "\t2\t"), "cache-hit field");
		badCache(cacheText.replace("g0001", "g0000"), "Duplicate cached genome");
		badCache(cacheText+"#CallGenesReplaySelectionCache\t1\n", "row width");
		models=load();
		cache=CallGenesReplaySelectionCache.readRelocatedCache(cacheFile, LABEL, models);
		write(mapping, mapText+"extra\n");
		try{
			models.verify();
			throw new AssertionError("Changed structural input accepted");
		}catch(IllegalStateException expected){
			check(expected.getMessage().contains("Structural input changed"), "Unexpected structural failure");
		}
		write(mapping, mapText);
		load().verify();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Mutation helpers       ----------------*/
	/*--------------------------------------------------------------*/

	/** Read current synthetic inventory/map with inherited labels and structural path checks. */
	private Inventory load()throws Exception{
		return CallGenesReplaySelectionCache.readRelocatedInventory(inventory, LABEL, mapping, resources);
	}
	/** Replace the synthetic map, require the named rejection, then restore its original text.
	 * Unexpected failures propagate before restoration so their input remains inspectable. */
	private void badMap(String text, String diagnostic)throws Exception{
		assert(text!=null) : "Each map mutation is an explicit fixture input";
		write(mapping, text);
		reject(()->load(), diagnostic);
		write(mapping, mapText);
	}
	/** Replace the synthetic cache and validate against a freshly loaded inventory.
	 * Restore the cache after the expected rejection; unexpected failures preserve the mutation. */
	private void badCache(String text, String diagnostic)throws Exception{
		assert(text!=null) : "Each cache mutation is isolated from accepted manifests";
		write(cacheFile, text);
		final Inventory models=load();
		reject(()->CallGenesReplaySelectionCache.readRelocatedCache(cacheFile, LABEL, models), diagnostic);
		write(cacheFile, cacheText);
	}
	/** Write or replace one owned synthetic UTF-8 file; never used for retained strict inputs. */
	private static void write(Path path, String text)throws Exception{
		assert(path!=null && text!=null) : "Fixture writes require explicit owned paths and content";
		Files.write(path, text.getBytes(StandardCharsets.UTF_8));
	}
	/** Change the first hexadecimal digit of a retained 20-character identity without hashing. */
	private static String wrongPin(String pin){
		check(pin.length()==20, "Existing control identity width");
		return (pin.charAt(0)=='0' ? "1" : "0")+pin.substring(1);
	}
	/** Require IllegalArgumentException containing the named diagnostic; count that check.
	 * A normal return, unrelated exception or mismatched diagnostic fails the fixture. */
	private static void reject(Action action, String diagnostic)throws Exception{
		assert(action!=null && diagnostic!=null) : "Negative fixtures require the intended failure, not an incidental I/O error";
		try{
			action.run();
		}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(diagnostic), "Unexpected rejection: "+expected.getMessage());
			return;
		}
		throw new AssertionError("Expected rejection: "+diagnostic);
	}
	/** Increment the process-local successful-check count or fail immediately. */
	private static void check(boolean pass, String message){
		if(!pass){throw new AssertionError(message);}
		checks++;
	}
	/** One potentially throwing fixture action used by the diagnostic-specific rejection helper. */
	private interface Action{
		/** Execute the negative case; unexpected checked failures propagate. */
		void run()throws Exception;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Fields and constants     ----------------*/
	/*--------------------------------------------------------------*/

	/** Owned synthetic resource directory, old/new payload paths and mutable metadata files. */
	private final Path resources, oldA, oldB, newA, newB, inventory, mapping, cacheFile;
	/** Synthetic FASTA name and immutable original metadata text used for restoration. */
	private final String fasta, inventoryText, mapText, cacheText;
	/** Synthetic inherited labels for fixture-only metadata; not measured digests. */
	private static final String LABEL="0123456789abcdef0123", OTHER="fedcba9876543210fedc";
	/** Exact version-one relocation magic and ordered column names. */
	private static final String MAP_HEADER="#CallGenesReplayModelRelocation\t1\noriginal_path\trelative_path\n";
	/** Process-local passed-check count; main() does not reset it for repeated invocation. */
	private static int checks;
}
