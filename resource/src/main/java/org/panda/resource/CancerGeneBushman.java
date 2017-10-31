package org.panda.resource;

import org.panda.utility.CollectionUtil;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * A curated cancer gene list from Bushman lab.
 * From http://www.bushmanlab.org/links/genelists
 *
 * @author Ozgun Babur
 */
public class CancerGeneBushman extends FileServer
{
	private static CancerGeneBushman instance;

	private Map<String, Set<String>> sym2resource;
	private Map<String, Set<String>> resource2sym;

	public static synchronized CancerGeneBushman get()
	{
		if (instance == null) instance = new CancerGeneBushman();
		return instance;
	}

	public Set<String> getAllSymbols()
	{
		return sym2resource.keySet();
	}

	public boolean isCancerGene(String sym)
	{
		return sym2resource.containsKey(sym);
	}

	public Set<String> getAllResources()
	{
		return resource2sym.keySet();
	}

	public Set<String> getSubset(String resource)
	{
		return resource2sym.get(resource);
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"cancer-genes-bushman.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://www.bushmanlab.org/assets/doc/allonco_20130923.tsv"};
	}

	@Override
	public boolean load() throws IOException
	{
		sym2resource = new HashMap<>();
		resource2sym = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0]).skip(1).map(l -> l.split("\t")).filter(t -> t[5].equals("human"))
			.forEach(token -> {
				String sym = token[2];
				String res = token[7];

				if (!sym2resource.containsKey(sym)) sym2resource.put(sym, new HashSet<>());
				sym2resource.get(sym).add(res);
				if (!resource2sym.containsKey(res)) resource2sym.put(res, new HashSet<>());
				resource2sym.get(res).add(sym);
			});

		return true;
	}

	public static void main(String[] args)
	{
		Set<String> okb = OncoKB.get().getAllSymbols();
		Set<String> bush = get().getAllSymbols();
		Set<String> cgc = CancerGeneCensus.get().getAllSymbols();
		CollectionUtil.printNameMapping("oncokb", "bushman", "cgc");
		CollectionUtil.printVennCounts(okb, bush, cgc);

		for (String res : get().getAllResources())
		{
			System.out.println();
			Set<String> sub = get().getSubset(res);
			CollectionUtil.printNameMapping(res, "cgc");
			CollectionUtil.printVennCounts(sub, cgc);
		}
	}
}
