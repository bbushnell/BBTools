package shared;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.HashMap;
import java.util.Map;

/**
 * Utility to automatically copy matching fields from a source to a dest.
 * Designed to reduce boilerplate in constructors.
 * 
 * @author Brian Bushnell
 * @contributor Sucrose
 * @date February 8, 2026
 */
public class FieldCopier {
	
	/** 
	 * Copies primitives and Strings from source to dest. 
	 * Ignores static fields and final fields in dest.
	 */
	public static int copy(Object dest, Object source) {
		return copy(dest, source, false, false, true);
	}

	/**
	 * Copies values of fields from source to dest where name and type match.
	 * @param dest The object to receive values
	 * @param source The object providing values
	 * @param includeObjects If true, copies complex Objects (Arrays, Lists, etc). 
	 * If false, only copies Primitives and Strings.
	 * @param includeStatics If true, also attempts to copy static fields.
	 * @return Number of fields set
	 */
	public static int copy(Object dest, Object source, 
			boolean includeObjects, boolean includeStatics, boolean includeSuperclass) {
		if(dest==null || source==null) {return 0;}
		
		// Map destination fields for O(1) lookup
		Map<String, Field> destFields=getAllFields(dest.getClass(), includeSuperclass);
		Class<?> sourceClass=source.getClass();
		int matches=0;
		
		while(sourceClass!=null && sourceClass!=Object.class) {
			Field[] fields=sourceClass.getDeclaredFields();
			
			for(Field srcField : fields) {
				try {
					String name=srcField.getName();
					int srcMods=srcField.getModifiers();
					
					// 1. Source Filter: Skip Statics unless requested
					if(Modifier.isStatic(srcMods) && !includeStatics) {continue;}
					
					// 2. Source Filter: Skip Complex Objects unless requested
					Class<?> type=srcField.getType();
					boolean isPrimOrString=(type.isPrimitive() || type==String.class);
					if(!isPrimOrString && !includeObjects) {continue;}
					
					// 3. Find Destination Field
					Field destField=destFields.get(name);
					if(destField==null) {continue;}
					
					// 4. Destination Filter: SAFETY CHECK - Skip Final Fields
					// We do not want to mess with JIT inlining or Definite Assignment rules.
					if(Modifier.isFinal(destField.getModifiers())) {continue;}

					// 5. Compatibility Check
					if(destField.getType().isAssignableFrom(type)) {
						
						// Handle accessibility
						if(!srcField.isAccessible()) {srcField.setAccessible(true);}
						if(!destField.isAccessible()) {destField.setAccessible(true);}
						
						// The Copy
						Object val=srcField.get(source);
						destField.set(dest, val);
						matches++;
					}
				} catch (Exception e) {
					// Security exceptions or access violations are ignored
				}
			}
			sourceClass=sourceClass.getSuperclass();
		}
		return matches;
	}
	
	/**
	 * Recursive helper to get all fields up the hierarchy.
	 * Handles shadowing by prioritizing the lowest class in the hierarchy.
	 */
	private static Map<String, Field> getAllFields(Class<?> clazz, boolean recursive){
		Map<String, Field> map=new HashMap<>();
		while(clazz!=null && clazz!=Object.class) {
			for(Field f : clazz.getDeclaredFields()) {
				// putIfAbsent ensures that if a subclass shadows a field, 
				// we keep the subclass version (which is usually what we want).
				map.putIfAbsent(f.getName(), f);
			}
			clazz=(recursive ? clazz.getSuperclass() : null);
		}
		return map;
	}
}