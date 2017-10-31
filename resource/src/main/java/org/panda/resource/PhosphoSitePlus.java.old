package org.panda.resource;

import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.TermCounter;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Provides phospho-site effects from PhosphoSitePlus.
 *
 * @author Ozgun Babur
 */
public class PhosphoSitePlus extends FileServer
{
	private static PhosphoSitePlus instance;

	Map<String, Map<String, Integer>> typeMap;
	Map<String, Map<String, String>> actualMap;

	public static synchronized PhosphoSitePlus get()
	{
		if (instance == null) instance = new PhosphoSitePlus();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"Regulatory_sites", "manually-curated-sites.txt"};
	}

	public Integer getEffect(String gene, String site)
	{
		if (typeMap.containsKey(gene))
		{
			return typeMap.get(gene).get(site);
		}
		return null;
	}

	public Integer getClosestEffect(String gene, String site, int distanceThreshold)
	{
		if (typeMap.containsKey(gene))
		{
			int s0 = Integer.parseInt(site.substring(1));

			Integer effect = null;
			int closestDist = Integer.MAX_VALUE;

			for (String ss : typeMap.get(gene).keySet())
			{
				Integer eff = typeMap.get(gene).get(ss);

				int s1 = Integer.parseInt(ss.substring(1));

				int dist = Math.abs(s0 - s1);
				if (dist == closestDist && !eff.equals(effect))
				{
					effect = 0;
				}
				else if (dist < closestDist)
				{
					effect = eff;
					closestDist = dist;
				}
			}
			if (closestDist <= distanceThreshold) return effect;
		}
		return null;
	}

	private List<String> sortSites(Set<String> sites)
	{
		List<String> list = new ArrayList<>(sites);
		Collections.sort(list, (o1, o2) -> {
			try
			{
				return new Integer(o1.substring(1)).compareTo(new Integer(o2.substring(1)));
			}
			catch (NumberFormatException e)
			{
				return 0;
			}
		});
		return list;
	}

	private List<String> getGenesWithMostSites()
	{
		List<String> genes = new ArrayList<>(typeMap.keySet());
		Collections.sort(genes, (o1, o2) -> new Integer(typeMap.get(o2).size()).compareTo(typeMap.get(o1).size()));
		return genes;
	}

	private void printSites(String gene)
	{
		System.out.println("Gene: " + gene);
		if (typeMap.containsKey(gene))
		{
			for (String site : sortSites(typeMap.get(gene).keySet()))
			{
				Integer sign = typeMap.get(gene).get(site);
				System.out.print("\tsite: " + site + "\t" + (sign == 1 ? "activating" : sign == -1 ? "inhibiting" : "complex"));
				System.out.println("\t(" + actualMap.get(gene).get(site) + ")");
			}
		}
	}

	@Override
	public boolean load() throws IOException
	{
		typeMap = new HashMap<>();
		actualMap = new HashMap<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[1]))).filter(l -> !l.startsWith("#"))
			.map(line -> line.split("\\s+")).filter(token -> token.length > 2).forEach(token ->
		{
			String gene = token[0];

			if (!typeMap.containsKey(gene)) typeMap.put(gene, new HashMap<>());
			if (!actualMap.containsKey(gene)) actualMap.put(gene, new HashMap<>());

			String site = token[1];
			int sign = Integer.parseInt(token[2]);

			typeMap.get(gene).put(site, sign);
			actualMap.get(gene).put(site, "manual curation");
		});

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0])), Charset.forName("windows-31j")).skip(4)
			.map(line -> line.split("\t"))
			.filter(token -> token.length >= 13 && token[6].equals("human") &&
				token[8].equals("PHOSPHORYLATION") && HGNC.get().getSymbol(token[4]) != null)
			.forEach(token -> {
				String gene = HGNC.get().getSymbol(token[4]);
				if (!typeMap.containsKey(gene)) typeMap.put(gene, new HashMap<>());
				if (!actualMap.containsKey(gene)) actualMap.put(gene, new HashMap<>());

				String site = token[7];

				if (typeMap.get(gene).containsKey(site)) return;

				actualMap.get(gene).put(site, token[12]);

				boolean actWord = false;
				boolean inhWord = false;

				if ((token[12].contains("induced") &&
					!token[12].contains("receptor desensitization, induced")))
				{
					actWord = true;
				}
				if ((token[12].contains("inhibited") ||
					token[12].contains("receptor desensitization, induced")))
				{
					inhWord = true;
				}

				if (actWord == inhWord)
				{
					if (token[12].contains("stabilization"))
					{
						actWord = true;
					}
					if (token[12].contains("degradation"))
					{
						inhWord = true;
					}
				}

				if (actWord == inhWord)
				{
					typeMap.get(gene).put(site, 0);
				} else
				{
					typeMap.get(gene).put(site, actWord ? 1 : -1);
				}
		});

		return true;
	}

	void printUniqueAA()
	{
		Set<String> sites = new HashSet<>();
		for (String gene : typeMap.keySet())
		{
			sites.addAll(typeMap.get(gene).keySet());
		}
		TermCounter tc = new TermCounter();
		for (String site : sites)
		{
			tc.addTerm(site.substring(0, 1));
		}
		tc.print();
	}

	public void fillInMissingEffect(Collection<ProteomicsFileRow> datas, int proximityThreshold)
	{
		for (ProteomicsFileRow data : datas)
		{
			if (data.effect != null) continue;
			if (data.sites == null || data.sites.isEmpty()) continue;

			Set<Integer> found = getEffects(data, proximityThreshold);
			data.effect = aggregateEffects(found);
		}
	}

	public Set<Integer> getEffects(ProteomicsFileRow data, int proximityThreshold)
	{
		Set<Integer> found = new HashSet<>();

		for (String gene : data.sites.keySet())
		{
			for (String site : data.sites.get(gene))
			{
				Integer e = getEffect(gene, site);
				if (e != null) found.add(e);
			}
		}

		if (found.isEmpty() && proximityThreshold > 0)
		{
			for (String gene : data.sites.keySet())
			{
				for (String site : data.sites.get(gene))
				{
					Integer e = getClosestEffect(gene, site, proximityThreshold);
					if (e != null) found.add(e);
				}
			}
		}

		return found;
	}

	private ProteomicsFileRow.SiteEffect aggregateEffects(Set<Integer> found)
	{
		if (found.contains(1))
		{
			if (found.contains(-1)) return ProteomicsFileRow.SiteEffect.COMPLEX;
			else return ProteomicsFileRow.SiteEffect.ACTIVATING;
		}
		else if (found.contains(-1))
		{
			return ProteomicsFileRow.SiteEffect.INHIBITING;
		}
		else if (!found.isEmpty()) return ProteomicsFileRow.SiteEffect.COMPLEX;
		return null;
	}

	public static void main(String[] args)
	{
		PhosphoSitePlus psp = new PhosphoSitePlus();
//		List<String> list = getGenesWithMostSites();
//		for (int i = 0; i < 10; i++)
//		{
//			printSites(list.get(i));
//		}
		psp.printSites("SPRY2");
//		printUniqueAA();

//		List<Integer> dists = new ArrayList<>();
//		for (String gene : typeMap.keySet())
//		{
//			Map<String, Integer> sites = typeMap.get(gene);
//			int min = Integer.MAX_VALUE;
//
//			for (String s1 : sites.keySet())
//			{
//				for (String s2 : sites.keySet())
//				{
//					if (sites.get(s1) * sites.get(s2) == -1)
//					{
//						int dif = Math.abs(Integer.parseInt(s1.substring(1)) - Integer.parseInt(s2.substring(1)));
//						if (dif < min) min = dif;
//					}
//				}
//			}
//
//			if (min < Integer.MAX_VALUE) dists.add(min);
//			if (min < 10)
//			{
//				System.out.println("\n" + gene + "\t" + min);
//				printSites(gene);
//			}
//		}
//
//		Histogram h = new Histogram(10);
//		h.setBorderAtZero(true);
//		for (Integer dist : dists)
//		{
//			h.count(dist);
//		}
//		h.print();
	}

}
