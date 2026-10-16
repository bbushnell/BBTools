package prok;

import java.io.File;
import dna.Data;

/** Sole owner of family-directory and member paths. No I/O occurs until resolve.
 * Loose files only; a future jar backend must supply reader capabilities explicitly.
 * @author Brian, Raiden */
public final class FamilyResourceResolver {

	public Handle resolve(String familyId, String role, Overrides overrides){
		final Handle h=resolve(familyId, role, ".", ".", overrides);
		if(role.equals("tables")){h.resolver=this;h.overrides=overrides;}return h;
	}
	private Handle resolve(String familyId, String role, String model, String end, Overrides overrides){
		final NcrnaFamilyRegistry.Entry family=NcrnaFamilyRegistry.get(familyId);
		final Overrides o=overrides==null ? new Overrides() : overrides;
		final String override=o.path(role);
		if(override!=null){return checked(family.id, role, override, role.equals("endpoints"));}
		if(role.equals("table") && o.tables!=null){return checked(family.id, role, new File(o.tables, model+'.'+end+".tsv.gz").toString(), false);}
		final String manifestPath=o.root==null ? Data.findPath("?"+family.directory+"/manifest.tsv", false) : new File(o.root, "manifest.tsv").toString();
		require(manifestPath!=null && new File(manifestPath).isFile(), "Missing required "+family.id+" manifest: "+manifestPath+" ("+family.directory+")");
		if(!manifestPath.equals(cachedPath)){
			manifest=new FamilyResourceManifest(manifestPath, family);cachedPath=manifestPath;
		}
		final FamilyResourceManifest.Entry entry=manifest.get(role, model, end);
		final Handle handle=checked(family.id, role, new File(new File(manifestPath).getParentFile(), entry.path).toString(), false);
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
		/** Root names one family directory, not the global resources directory. */
		public String root, consensus, models, seeds, tables, endpoints;
		private String path(String role){
			if(role.equals("consensus")){return consensus;}if(role.equals("models")){return models;}
			if(role.equals("seeds")){return seeds;}if(role.equals("endpoints")){return endpoints;}
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
	private String cachedPath;
	private FamilyResourceManifest manifest;
}
