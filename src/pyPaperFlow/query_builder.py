import typer
import requests
import json
from typing import List, Dict, Optional
from datetime import datetime

class AIQueryAssistant:
    def __init__(self, api_key: str, api_base: str = "https://api.openai.com/v1", model: str = "gpt-4o"):
        self.api_key = api_key
        self.api_base = api_base
        self.model = model

    def chat(self, messages: List[Dict[str, str]]) -> str:
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }
        data = {
            "model": self.model,
            "messages": messages,
            "temperature": 0.7
        }
        try:
            response = requests.post(f"{self.api_base}/chat/completions", headers=headers, json=data)
            response.raise_for_status()
            return response.json()["choices"][0]["message"]["content"]
        except Exception as e:
            typer.echo(f"AI API Error: {e}")
            return ""

    def refine_query_plan(self, current_info: Dict[str, str]) -> Dict[str, any]:
        """
        Ask AI to analyze current inputs, suggest MeSH terms, synonyms, and next questions.
        """
        system_prompt = """You are a PubMed Search Expert. Your goal is to help the user construct a high-quality, precise, and comprehensive PubMed search query.
        
        You have the following tasks:
        1. Analyze the user's current inputs.
        2. Identify relevant MeSH terms and keywords (including synonyms, spelling variations).
        3. Suggest logical operators (AND, OR, NOT) to combine these terms.
        4. Identify missing information or potential ambiguities.
        5. Propose 1-3 follow-up questions to refine the search scope.
        
        Output your response in JSON format with the following structure:
        {
            "analysis": "Brief analysis of the current search intent.",
            "suggestions": {
                "mesh_terms": ["Term1", "Term2"],
                "keywords": ["Keyword1", "Keyword2"],
                "logic_advice": "Advice on how to combine terms."
            },
            "draft_query": "A draft PubMed query string based on current info.",
            "follow_up_questions": ["Question 1", "Question 2"]
        }
        """
        
        user_content = f"Current User Inputs: {json.dumps(current_info)}"
        
        response_text = self.chat([
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_content}
        ])
        
        try:
            # Try to parse JSON from the response (handle potential markdown code blocks)
            if "```json" in response_text:
                response_text = response_text.split("```json")[1].split("```")[0].strip()
            elif "```" in response_text:
                response_text = response_text.split("```")[1].split("```")[0].strip()
            return json.loads(response_text)
        except json.JSONDecodeError:
            return {"error": "Failed to parse AI response", "raw_response": response_text}

class QueryBuilder:
    def __init__(self, ai_assistant: Optional[AIQueryAssistant] = None):
        self.ai = ai_assistant
        self.inputs = {}
        self.history = []

    def ask(self, question: str, key: str, default: str = ""):
        answer = typer.prompt(question, default=default)
        if answer:
            self.inputs[key] = answer
        return answer

    def run(self):
        typer.secho("Welcome to the PubMed Query Builder Wizard!", fg=typer.colors.GREEN, bold=True)
        typer.echo("I will guide you through constructing a complex search query step-by-step.")
        
        # Step 1: Core Concepts
        typer.secho("\n--- Step 1: Core Concepts ---", fg=typer.colors.CYAN)
        self.ask("What is the main subject/gene/protein/entity you are interested in? (e.g., CTCF, p53)", "subject")
        self.ask("Are there any specific diseases, conditions, or biological processes involved? (e.g., Cancer, Folding)", "process_condition")
        self.ask("Are you looking for a specific method or technique? (e.g., ChIP-seq, AlphaFold)", "method")
        
        # AI Refinement 1
        if self.ai:
            typer.secho("\n... Consulting AI Expert for refinement ...", fg=typer.colors.MAGENTA)
            advice = self.ai.refine_query_plan(self.inputs)
            if "error" not in advice:
                self.display_ai_advice(advice)
                
                # Ask follow-up questions suggested by AI
                if "follow_up_questions" in advice and advice["follow_up_questions"]:
                    typer.secho("\n--- AI Suggested Follow-up ---", fg=typer.colors.CYAN)
                    for i, q in enumerate(advice["follow_up_questions"]):
                        ans = typer.prompt(f"AI Q{i+1}: {q}", default="")
                        if ans:
                            self.inputs[f"ai_q_{i}"] = f"Q: {q} A: {ans}"
            else:
                typer.echo(f"AI Error: {advice.get('raw_response')}")

        # Step 2: Constraints & Exclusions
        typer.secho("\n--- Step 2: Constraints & Exclusions ---", fg=typer.colors.CYAN)
        self.ask("Any specific organism? (e.g., Human, Mouse)", "organism")
        self.ask("Any keywords to EXCLUDE? (e.g., Review, specific authors)", "exclusions")
        self.ask("Date range? (e.g., 2020:2025, last 5 years)", "date_range")
        
        # Final Construction
        typer.secho("\n--- Finalizing Query ---", fg=typer.colors.CYAN)
        final_query = self.construct_query()
        
        if self.ai:
            typer.secho("\n... Asking AI to polish the final query ...", fg=typer.colors.MAGENTA)
            final_polish = self.ai.chat([
                {"role": "system", "content": "You are a PubMed Search Expert. Refine the following PubMed query to be syntactically correct and optimized. Output ONLY the query string."},
                {"role": "user", "content": f"User Inputs: {json.dumps(self.inputs)}\nDraft Query: {final_query}"}
            ])
            if final_polish:
                final_query = final_polish.strip().strip('"').strip("'")

        typer.secho("\nGenerated PubMed Query:", fg=typer.colors.GREEN, bold=True)
        typer.echo(final_query)
        
        return final_query

    def display_ai_advice(self, advice):
        typer.secho("\n[AI Analysis]", fg=typer.colors.YELLOW)
        typer.echo(advice.get("analysis", ""))
        
        suggestions = advice.get("suggestions", {})
        if suggestions.get("mesh_terms"):
            typer.echo(f"Suggested MeSH Terms: {', '.join(suggestions['mesh_terms'])}")
        if suggestions.get("keywords"):
            typer.echo(f"Suggested Keywords: {', '.join(suggestions['keywords'])}")
        if advice.get("draft_query"):
            typer.echo(f"Draft Query Idea: {advice['draft_query']}")

    def construct_query(self):
        # Basic rule-based construction if AI fails or as a base
        parts = []
        
        if self.inputs.get("subject"):
            parts.append(f"({self.inputs['subject']}[Title/Abstract])")
            
        if self.inputs.get("process_condition"):
            parts.append(f"({self.inputs['process_condition']}[Title/Abstract])")
            
        if self.inputs.get("method"):
            parts.append(f"({self.inputs['method']}[Title/Abstract])")
            
        if self.inputs.get("organism"):
            parts.append(f"({self.inputs['organism']}[Title/Abstract])")
            
        base_query = " AND ".join(parts)
        
        if self.inputs.get("exclusions"):
            base_query += f" NOT ({self.inputs['exclusions']})"
            
        if self.inputs.get("date_range"):
            # Simple date handling, user might need to format it correctly or AI handles it
            base_query += f" AND ({self.inputs['date_range']}[Date - Publication])"
            
        return base_query
