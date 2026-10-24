package prok;

import java.io.File;
import dna.Data;

/** Sole owner of family-directory and member paths. No I/O occurs until resolve.
 * Loose files only; a future jar backend must supply reader capabilities explicitly.
 * @author Brian, Raiden */
public final class FamilyResourceResolver {

	/** Retain a declared private local reference for legacy tools; published installs
	 * use the canonical family consensus. The production family loader never calls this. */
	public Handle legacyConsensus(String familyId){
		final NcrnaFamilyRegistry.Entry family=NcrnaFamilyRegistry.get(familyId);
		if(family.privateLegacyConsensus!=null){
			final File local=new File(Data.RESOURCES(), family.privateLegacyConsensus);
			if(local.isFile()){return checked(family.id, "consensus", local.toString(), false);}
		}return resolve(family.id, family.legacyConsensusRole, null);
	}
	/** Explicit historical tRNA K values retain their flat compatibility resources. */
	Handle seeds(String familyId, int k, Overrides overrides){
		final NcrnaFamilyRegistry.Entry family=NcrnaFamilyRegistry.get(familyId);
		if(k==family.seedK){return resolve(family.id, "seeds", overrides);}
		require(k>0 && k<=31 && family.boundary!=null && family.boundary.trna,
			"Only the independent tRNA caller supports explicit legacy seed K: "+familyId+"/"+k);
		if(overrides!=null && overrides.seeds!=null){return checked(family.id, "seeds", overrides.seeds, false);}
		final String stem=family.id+'_'+k+"mers.fa";
		String path=Data.findPath("?"+stem, false);
		if(path==null){path=Data.findPath("?"+stem+".gz", false);}
		require(path!=null, "Missing explicit legacy "+family.id+" seeds: "+stem);
		return checked(family.id, "seeds", path, false);
	}

	public Handle resolve(String familyId, String role, Overrides overrides){
		final Handle h=resolve(familyId, role, ".", ".", overrides);
		if(role.equals("tables")){h.resolver=this;h.overrides=overrides;}return h;
	}
	private Handle resolve(String familyId, String role, String model, String end, Overrides overrides){
		final NcrnaFamilyRegistry.Entry family=NcrnaFamilyRegistry.get(familyId);
		final Overrides o=overrides==null ? new Overrides() : overrides;
		final String override=o.path(role);
		if(override!=null){
			// Modern18S historically accepted literal paths and resource-relative names.
			final String chosen=family.resources!=null && family.resources.cwdName(role)!=null
				&& !new File(override).isFile() && !override.startsWith("?") ? "?"+override : override;
			return checked(family.id, role, chosen, role.equals("endpoints"));
		}
		if(o.root==null && family.resources!=null){
			final String legacy=family.resources.cwdName(role);
			if(legacy!=null && new File(legacy).isFile()){return checked(family.id, role, legacy, false);}
		}
		// The old generic factories explicitly preferred a literal CWD consensus/plain HBM.
		// Explicit roots and overrides bypass this compatibility lookup entirely.
		if(o.root==null && family.boundary!=null){
			final String legacy=role.equals("consensus") ? family.boundary.consensus : role.equals("models") ? family.boundary.models : null;
			if(legacy!=null && new File(legacy).isFile()){return checked(family.id, role, legacy, false);}
		}
		if(role.equals("table") && o.tables!=null){return checked(family.id, role, new File(o.tables, model+'.'+end+".tsv.gz").toString(), false);}
		final String manifestPath=o.root==null ? Data.findPath("?"+family.directory+"/manifest.tsv", false) : new File(o.root, "manifest.tsv").toString();
		require(manifestPath!=null && new File(manifestPath).isFile(), "Missing required "+family.id+" manifest: "+manifestPath+" ("+family.directory+")");
		final String identity=family.group==null ? "1:"+family.id : "2:"+family.group.id;
		// Members share their group's cache; different families must still validate the header at an explicit root.
		if(!manifestPath.equals(cachedPath) || !identity.equals(cachedIdentity)){
			manifest=new FamilyResourceManifest(manifestPath, family);cachedPath=manifestPath;cachedIdentity=identity;
		}
		final String member=family.group!=null && family.group.shared(role) ? "." : family.member;
		final FamilyResourceManifest.Entry entry=manifest.get(member, role, model, end);
		File directory=new File(manifestPath).getParentFile();
		if(!entry.member.isEmpty() && !entry.member.equals(".")){directory=new File(directory, entry.member);}
		final Handle handle=checked(family.id, role, new File(directory, entry.path).toString(), false);
		require(new File(handle.path).length()==entry.bytes, "Incomplete "+family.id+" resource "+role+": "+handle.path+"; expected "+entry.bytes+" bytes");
		return handle;
	}
	public Handle tables(String familyId, Overrides overrides){
		return resolve(familyId, "tables", overrides);
	}
	/** Explicit resources keep their paths; no fallback to a different installation. */
	private static Handle checked(String family, String role, String path, boolean directoryAllowed){
		final String resolved=path.startsWith("?") ? Data.findPath(path, false) : path;
		require(resolved!=null && (new File(resolved).isFile() || directoryAllowed && new File(resolved).isDirectory()),
			"Missing required "+family+" resource "+role+": "+path);
		return new Handle(family, role, resolved);
	}
	public static final class Overrides {
		/** Root names one family directory (the group root for grouped members), not the global resources directory. */
		public String root, consensus, models, seeds, tables, endpoints, cutoffs;
		public String boundary, boundary3, boundaryStart, boundaryStop;
		public String voteOffsets, joinedModelPositions, lengthBounds;
		private String path(String role){
			if(role.equals("consensus")){return consensus;}if(role.equals("models")){return models;}
			if(role.equals("seeds")){return seeds;}if(role.equals("endpoints")){return endpoints;}
			if(role.equals("cutoffs")){return cutoffs;}
			if(role.equals("length_bounds")){return lengthBounds;}
			if(role.equals("vote_offsets")){return voteOffsets;}
			if(role.equals("joined_model_positions")){return joinedModelPositions;}
			if(role.equals("boundary")){return boundary;}
			if(role.equals("boundary_start")){return boundaryStart;}
			if(role.equals("boundary_stop")){return boundaryStop;}
			if(role.equals("tables") && tables!=null){return new File(tables, "TABLES.tsv").toString();}return null;
		}
	}
	public static final class Handle {
		private Handle(String family_, String role_, String path_){family=family_;role=role_;path=path_;}
		public final String family, role, path;
		public boolean isFile(){return new File(path).isFile();}
		String table(String model, String end){
			assert(resolver!=null && role.equals("tables")) : "Table children must resolve through their family manifest";
			return resolver.resolve(family, "table", model, end, overrides).path;
		}
		String child(String relative){
			FamilyResourceManifest.validateRelative(relative);
			return checked(family, role, new File(path, relative).toString(), false).path;
		}
		String optionalChild(String relative){
			FamilyResourceManifest.validateRelative(relative);
			final String value=new File(path, relative).toString();return new File(value).isFile() ? value : null;
		}
		private FamilyResourceResolver resolver;
		private Overrides overrides;
	}
	private static void require(boolean ok, String why){FamilyResourceManifest.require(ok, why);}
	private String cachedPath, cachedIdentity;
	private FamilyResourceManifest manifest;
}
