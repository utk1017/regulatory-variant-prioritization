import argparse

def main():
    parser = argparse.ArgumentParser(description="Regulatory Variant Prioritization Framework")
    parser.add_argument("command", choices=["run"], help="Command to execute")
    parser.add_argument("--layer", type=int, required=True, help="Layer to run (e.g., 0)")
    parser.add_argument("--config", type=str, required=True, help="Path to config YAML")

    args = parser.parse_args()

    if args.command == "run":
        if args.layer == 0:
            print(f"Running Layer 0 with config {args.config}...")
            # from regvar.layer0.pipeline import run_layer0
            # run_layer0(args.config)
        else:
            print(f"Layer {args.layer} is not yet implemented.")

if __name__ == "__main__":
    main()
